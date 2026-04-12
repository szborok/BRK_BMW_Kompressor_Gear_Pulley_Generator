# =============================================================================
# Fusion 360 Python API Script
# BMW Kompressor gear pulley generator
# =============================================================================

import math
import os
import traceback
from bisect import bisect_left

import adsk.core
import adsk.fusion


# Reference DXF — must sit in the same folder as this script
REFERENCE_DXF = "Sketch14.dxf"
REFERENCE_TOOTH_COUNT = 57      # periodicity measured in Sketch14

# -----------------------------------------------------------------------
# Output gear — change only PROFILE_SCALE and TOOTH_DEPTH_SCALE to tune
# -----------------------------------------------------------------------
DEFAULT_TOOTH_COUNT = 44
FACE_WIDTH    = 30.0            # mm, extrusion depth

# >1.0 = deeper valleys / more prominent teeth  (was 1.0 = too shallow)
TOOTH_DEPTH_SCALE = 1.30

# Uniform XY scale — increase to make the whole gear bigger
PROFILE_SCALE = 0.7984

# Small pulley compensation:
# 0.0 = legacy uniform scaling (shallower teeth at low tooth counts)
# 1.0 = keep tooth depth offsets constant as tooth count changes
SMALL_PULLEY_DEPTH_COMP = 1.0

# Rotate so one tooth tip points to +Y
ALIGN_TOOTH_TO_POSITIVE_Y = True

VERSION = "4.1"

POINTS_PER_TOOTH = 20
MAX_PROFILE_POINTS = 900

# Small pulley pitch compression: <1.0 = tighter tooth spacing on low tooth counts.
SMALL_PULLEY_PITCH_COMP = 0.945

# Small pulley tooth width shaping: <1.0 = wider teeth / narrower gaps.
SMALL_PULLEY_TOOTH_WIDTH_EXP = 0.85


def mm_to_cm(v):
    return v * 0.1


# -----------------------------------------------------------------------
# DXF loading
# -----------------------------------------------------------------------
def _read_dxf_pairs(path):
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        raw = [ln.rstrip("\r\n") for ln in f]

    pairs = []
    i = 0
    while i + 1 < len(raw):
        pairs.append((raw[i].strip(), raw[i + 1].strip()))
        i += 2
    return pairs


def _extract_first_lwpolyline_vertices(path):
    pairs = _read_dxf_pairs(path)

    in_entities = False
    entity_pairs = []
    i = 0
    while i < len(pairs):
        code, value = pairs[i]

        if not in_entities:
            if code == "2" and value == "ENTITIES":
                in_entities = True
            i += 1
            continue

        if code == "0" and value == "ENDSEC":
            break

        if code == "0" and value == "LWPOLYLINE":
            entity_pairs = []
            i += 1
            while i < len(pairs):
                c, v = pairs[i]
                if c == "0":
                    break
                entity_pairs.append((c, v))
                i += 1
            break

        i += 1

    if not entity_pairs:
        raise RuntimeError("No LWPOLYLINE found in reference DXF")

    closed = False
    vertices = []
    cur = None

    for code, value in entity_pairs:
        if code == "70":
            flags = int(float(value))
            closed = (flags & 1) != 0
        elif code == "10":
            if cur is not None and cur["x"] is not None and cur["y"] is not None:
                vertices.append(cur)
            cur = {"x": float(value), "y": None, "bulge": 0.0}
        elif code == "20" and cur is not None:
            cur["y"] = float(value)
        elif code == "42" and cur is not None:
            cur["bulge"] = float(value)

    if cur is not None and cur["x"] is not None and cur["y"] is not None:
        vertices.append(cur)

    if len(vertices) < 2:
        raise RuntimeError("Reference LWPOLYLINE has insufficient vertices")

    return vertices, closed


# -----------------------------------------------------------------------
# Profile processing helpers
# -----------------------------------------------------------------------
def _densify_polyline(vertices_mm, closed, arc_steps=16):
    """Expand bulge arcs to dense (x,y) point list."""
    pts = []
    n = len(vertices_mm)
    seg_count = n if closed else n - 1
    for i in range(seg_count):
        v1 = vertices_mm[i]
        v2 = vertices_mm[(i + 1) % n]
        x1, y1 = v1["x"], v1["y"]
        x2, y2 = v2["x"], v2["y"]
        bulge = v1.get("bulge", 0.0)
        if i == 0:
            pts.append((x1, y1))
        if abs(bulge) < 1e-10:
            pts.append((x2, y2))
            continue
        theta = 4.0 * math.atan(bulge)
        dx, dy = x2 - x1, y2 - y1
        chord = math.hypot(dx, dy)
        if chord < 1e-12 or abs(math.tan(theta * 0.5)) < 1e-12:
            pts.append((x2, y2))
            continue
        mx, my = 0.5*(x1+x2), 0.5*(y1+y2)
        nx, ny = -dy/chord, dx/chord
        d  = chord / (2.0 * math.tan(theta * 0.5))
        cx, cy = mx + nx*d, my + ny*d
        rr = math.hypot(x1-cx, y1-cy)
        a0 = math.atan2(y1-cy, x1-cx)
        for j in range(1, arc_steps):
            t = j / float(arc_steps)
            a = a0 + theta * t
            pts.append((cx + rr*math.cos(a), cy + rr*math.sin(a)))
        pts.append((x2, y2))
    if closed and pts and math.hypot(pts[0][0]-pts[-1][0], pts[0][1]-pts[-1][1]) < 1e-9:
        pts.pop()
    return pts


def _polar_samples(pts):
    """Sorted (angle 0..2π, radius) for a ring of 2-D points."""
    cx = sum(p[0] for p in pts) / len(pts)
    cy = sum(p[1] for p in pts) / len(pts)
    arr = []
    for x, y in pts:
        a = math.atan2(y-cy, x-cx)
        if a < 0:
            a += 2.0 * math.pi
        arr.append((a, math.hypot(x-cx, y-cy)))
    arr.sort()
    return arr


def _interp_r(angles, radii, a):
    two_pi = 2.0 * math.pi
    a = a % two_pi
    n = len(angles)
    i = bisect_left(angles, a)
    if i == 0:
        a0, r0 = angles[-1] - two_pi, radii[-1]
        a1, r1 = angles[0],           radii[0]
    elif i >= n:
        a0, r0 = angles[-1],          radii[-1]
        a1, r1 = angles[0] + two_pi,  radii[0]
    else:
        a0, r0 = angles[i-1], radii[i-1]
        a1, r1 = angles[i],   radii[i]
    if abs(a1 - a0) < 1e-12:
        return r0
    return r0 + (r1 - r0) * (a - a0) / (a1 - a0)


def _find_first_tip_angle(angles, radii):
    """Angle of the tooth-tip center nearest to +X (angle 0)."""
    rmax  = max(radii)
    rmean = sum(radii) / len(radii)
    eps   = max(0.05, (rmax - rmean) * 0.15)
    tip   = sorted([(a, r) for a, r in zip(angles, radii) if r >= rmax - eps])
    if not tip:
        return max(zip(angles, radii), key=lambda t: t[1])[0]
    # cluster and find cluster closest to angle 0
    gap      = 2.0 * math.pi / len(angles) * 8
    clusters = [[tip[0][0]]]
    for a, _ in tip[1:]:
        if a - clusters[-1][-1] > gap:
            clusters.append([a])
        else:
            clusters[-1].append(a)
    best, best_d = clusters[0], 1e9
    for cl in clusters:
        m = sum(cl) / len(cl)
        d = abs(((m + math.pi) % (2*math.pi)) - math.pi)  # dist to 0
        if d < best_d:
            best_d, best = d, cl
    return sum(best) / len(best)


def _remap_teeth(polar, ref_teeth, new_teeth, depth_scale, out_pts, tooth_width_exp=1.0):
    """
    Re-tile a ref_teeth-periodic radial profile into new_teeth teeth.
    Amplifies radial oscillation by depth_scale (>1 = deeper/prominent).
    Returns list of (x,y) mm coordinates.
    """
    angles = [a for a, _ in polar]
    radii  = [r for _, r in polar]
    rmean  = sum(radii) / len(radii)

    T_ref = 2.0 * math.pi / ref_teeth

    tip0 = _find_first_tip_angle(angles, radii)   # center of tooth-0 in ref

    result = []
    for i in range(out_pts):
        a_out = (2.0 * math.pi * i) / out_pts

        # phase within current output tooth, range [-0.5 .. +0.5]
        phase_frac = math.fmod(a_out * new_teeth / (2.0 * math.pi), 1.0)
        if phase_frac > 0.5:
            phase_frac -= 1.0

        # Shape tooth width: exp < 1 widens teeth and narrows gaps.
        phase_norm = abs(phase_frac) / 0.5
        pow_norm = phase_norm ** max(1e-6, tooth_width_exp)

        local_ref = math.copysign(pow_norm * 0.5 * T_ref, phase_frac)

        # reference angle: tooth-0 center + local offset
        a_ref = tip0 + local_ref

        r_ref   = _interp_r(angles, radii, a_ref)
        r_final = rmean + (r_ref - rmean) * depth_scale

        result.append((r_final * math.cos(a_out), r_final * math.sin(a_out)))

    return result


def largest_profile(sketch):
    count = sketch.profiles.count
    if count == 0:
        return None
    if count == 1:
        return sketch.profiles.item(0)

    # Avoid expensive areaProperties() on complex sketches; use bounding box area.
    best = sketch.profiles.item(0)
    best_area = -1.0
    for i in range(count):
        pr = sketch.profiles.item(i)
        try:
            bb = pr.boundingBox
            area = (bb.maxPoint.x - bb.minPoint.x) * (bb.maxPoint.y - bb.minPoint.y)
        except Exception:
            area = -1.0
        if area > best_area:
            best_area = area
            best = pr
    return best


def _build_pts_for_tooth_count(polar, tooth_count):
    # Cap profile points to keep sketch generation responsive in Fusion.
    out_n = min(MAX_PROFILE_POINTS, tooth_count * POINTS_PER_TOOTH)

    if tooth_count < DEFAULT_TOOTH_COUNT:
        blend = (DEFAULT_TOOTH_COUNT - tooth_count) / float(DEFAULT_TOOTH_COUNT - 3)
        blend = max(0.0, min(1.0, blend))
    else:
        blend = 0.0

    width_exp = 1.0 + (SMALL_PULLEY_TOOTH_WIDTH_EXP - 1.0) * blend
    pts = _remap_teeth(
        polar,
        REFERENCE_TOOTH_COUNT,
        tooth_count,
        TOOTH_DEPTH_SCALE,
        out_n,
        width_exp,
    )
    if ALIGN_TOOTH_TO_POSITIVE_Y:
        cs, sn = math.cos(math.pi * 0.5), math.sin(math.pi * 0.5)
        pts = [(x*cs - y*sn, x*sn + y*cs) for x, y in pts]

    # Scale pitch radius with tooth count, and blend how much tooth depth scales.
    # This helps small pulleys seat the belt correctly while keeping 44T unchanged.
    tooth_scale = tooth_count / float(DEFAULT_TOOTH_COUNT)
    comp = max(0.0, min(1.0, SMALL_PULLEY_DEPTH_COMP))

    radii = [math.hypot(x, y) for x, y in pts]
    mean_radius = sum(radii) / len(radii)

    # Blend-in pitch compression only for smaller pulleys to tighten spacing.
    if tooth_count < DEFAULT_TOOTH_COUNT:
        pitch_scale = 1.0 - (1.0 - SMALL_PULLEY_PITCH_COMP) * blend
    else:
        pitch_scale = 1.0

    target_mean_radius = mean_radius * tooth_scale * pitch_scale

    # 0 -> legacy (depth scales with tooth_scale), 1 -> constant depth offsets.
    depth_scale = tooth_scale + comp * (1.0 - tooth_scale)

    scaled_pts = []
    for (x, y), radius in zip(pts, radii):
        if radius < 1e-12:
            scaled_pts.append((x, y))
            continue
        target_radius = target_mean_radius + (radius - mean_radius) * depth_scale
        scale = target_radius / radius
        scaled_pts.append((x * scale, y * scale))
    return scaled_pts


def _outer_diameter_for_tooth_count(polar, tooth_count, profile_scale):
    pts = _build_pts_for_tooth_count(polar, tooth_count)
    r_vals = [math.hypot(x, y) for x, y in pts]
    return 2 * max(r_vals) * profile_scale


# -----------------------------------------------------------------------
# Entry point
# -----------------------------------------------------------------------
def run(context):
    ui = None
    try:
        if adsk is None:
            raise RuntimeError("Run this script from Fusion 360 Scripts and Add-Ins")

        app = adsk.core.Application.get()
        ui  = app.userInterface
        design = adsk.fusion.Design.cast(app.activeProduct)
        root   = design.rootComponent
        title = f"BMW Kompressor gear pulley generator v{VERSION}"

        mode_text, mode_cancelled = ui.inputBox(
            "Choose mode: T = tooth count, D = desired diameter.\n"
            "Tip: you can also type a number directly (e.g. 33) to generate that tooth count.",
            title,
            "D"
        )
        if mode_cancelled:
            return
        mode_raw = mode_text.strip()
        mode = mode_raw.upper()

        direct_tooth_count = None
        if mode not in ("D", "T"):
            try:
                direct_tooth_count = int(mode_raw)
            except Exception:
                raise RuntimeError("Mode must be D (diameter), T (tooth count), or a direct integer tooth count")
            if direct_tooth_count < 3:
                raise RuntimeError("Tooth count must be >= 3")
            mode = "T"

        # ---- load & densify reference ----
        dxf_path = os.path.join(os.path.dirname(__file__), REFERENCE_DXF)
        if not os.path.exists(dxf_path):
            raise RuntimeError(f"Reference DXF not found: {dxf_path}")

        vertices, closed = _extract_first_lwpolyline_vertices(dxf_path)
        dense = _densify_polyline(vertices, closed, arc_steps=16)
        polar = _polar_samples(dense)

        profile_scale = PROFILE_SCALE

        if mode == "T":
            if direct_tooth_count is not None:
                tooth_count = direct_tooth_count
            else:
                tooth_text, tooth_cancelled = ui.inputBox(
                    "Tooth count:",
                    title,
                    str(DEFAULT_TOOTH_COUNT)
                )
                if tooth_cancelled:
                    return
                try:
                    tooth_count = int(tooth_text.strip())
                except Exception:
                    raise RuntimeError("Tooth count must be an integer")
                if tooth_count < 3:
                    raise RuntimeError("Tooth count must be >= 3")
        else:
            od_text, od_cancelled = ui.inputBox(
                "Desired outer diameter in mm:",
                title,
                "132"
            )
            if od_cancelled:
                return
            try:
                desired_outer_diameter = float(od_text.strip())
            except Exception:
                raise RuntimeError("Desired outer diameter must be a valid number")
            if desired_outer_diameter <= 0:
                raise RuntimeError("Desired outer diameter must be > 0")

            od_cache = {}

            def od_for(tc):
                if tc not in od_cache:
                    od_cache[tc] = _outer_diameter_for_tooth_count(polar, tc, profile_scale)
                return od_cache[tc]

            est_teeth = desired_outer_diameter / od_for(DEFAULT_TOOTH_COUNT) * DEFAULT_TOOTH_COUNT
            low_tc = max(3, int(math.floor(est_teeth)))
            high_tc = max(3, int(math.ceil(est_teeth)))
            if high_tc == low_tc:
                high_tc += 1

            low_od = od_for(low_tc)
            high_od = od_for(high_tc)

            # Expand outward until low <= desired <= high
            while low_tc > 3 and low_od > desired_outer_diameter:
                high_tc, high_od = low_tc, low_od
                low_tc -= 1
                low_od = od_for(low_tc)

            while high_od < desired_outer_diameter and high_tc < 300:
                low_tc, low_od = high_tc, high_od
                high_tc += 1
                high_od = od_for(high_tc)

            choose_text, choose_cancelled = ui.inputBox(
                "Closest options to desired diameter:\n\n"
                f"Desired OD: {desired_outer_diameter:.3f} mm\n"
                f"1) Smaller: {low_tc} teeth -> {low_od:.3f} mm\n"
                f"2) Bigger : {high_tc} teeth -> {high_od:.3f} mm\n\n"
                "Type 1 or 2",
                title,
                "1"
            )
            if choose_cancelled:
                return
            choice = choose_text.strip()
            if choice == "1":
                tooth_count = low_tc
            elif choice == "2":
                tooth_count = high_tc
            else:
                raise RuntimeError("Selection must be 1 or 2")

        # ---- remap to requested teeth with depth amplification ----
        pts = _build_pts_for_tooth_count(polar, tooth_count)

        raw_r_vals = [math.hypot(x, y) for x, y in pts]
        raw_tip_diameter = 2 * max(raw_r_vals)
        raw_root_diameter = 2 * min(raw_r_vals)

        preview_tip = raw_tip_diameter * profile_scale
        preview_root = raw_root_diameter * profile_scale
        preview_depth = preview_tip - preview_root

        preview_msg = (
            "Preview before generation:\n\n"
            f"Teeth: {tooth_count}\n"
            f"Outer diameter: {preview_tip:.3f} mm\n"
            f"Root diameter : {preview_root:.3f} mm\n"
            f"Tooth depth   : {preview_depth:.3f} mm\n"
            f"Depth scale   : {TOOTH_DEPTH_SCALE}\n"
            f"Profile scale : {profile_scale:.6f}\n"
            f"Face width    : {FACE_WIDTH:.3f} mm\n\n"
            "Type YES to generate the pulley, anything else cancels."
        )
        confirm_text, confirm_cancelled = ui.inputBox(preview_msg, title, "YES")
        if confirm_cancelled or confirm_text.strip().upper() != "YES":
            return

        # ---- apply selected scale ----
        pts = [(x * profile_scale, y * profile_scale) for x, y in pts]

        # ---- build sketch as smooth spline contour ----
        sketch = root.sketches.add(root.xYConstructionPlane)
        sketch.name = f"BMW_Kompressor_{tooth_count}T_v{VERSION}"
        sketch_pts = [adsk.core.Point3D.create(mm_to_cm(x), mm_to_cm(y), 0) for x, y in pts]

        try:
            fit_pts = adsk.core.ObjectCollection.create()
            for p in sketch_pts:
                fit_pts.add(p)
            # Repeat first point to enforce closure for the spline.
            fit_pts.add(sketch_pts[0])
            spline = sketch.sketchCurves.sketchFittedSplines.add(fit_pts)
            try:
                spline.isClosed = True
            except Exception:
                pass
        except Exception:
            # Fallback: draw as polyline if spline creation fails.
            lines = sketch.sketchCurves.sketchLines
            for i in range(len(sketch_pts)):
                p1 = sketch_pts[i]
                p2 = sketch_pts[(i + 1) % len(sketch_pts)]
                lines.addByTwoPoints(p1, p2)

        pr = largest_profile(sketch)
        if pr is None:
            raise RuntimeError("No closed profile found")

        extrudes = root.features.extrudeFeatures
        ext_in = extrudes.createInput(pr, adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
        ext_in.setDistanceExtent(False, adsk.core.ValueInput.createByReal(mm_to_cm(FACE_WIDTH)))
        body = extrudes.add(ext_in).bodies.item(0)
        body.name = f"BMW_Kompressor_{tooth_count}T_v{VERSION}"

        r_vals = [math.hypot(x, y) for x, y in pts]
        ui.messageBox(
            "Pulley generated.\n\n"
            f"Teeth: {tooth_count}\n"
            f"Tip  diameter : {2*max(r_vals):.3f} mm\n"
            f"Root diameter : {2*min(r_vals):.3f} mm\n"
            f"Tooth depth   : {max(r_vals)-min(r_vals):.3f} mm\n"
            f"Depth scale   : {TOOTH_DEPTH_SCALE}\n"
            f"Profile scale : {profile_scale:.6f}\n"
            f"Face width    : {FACE_WIDTH:.3f} mm"
        )

    except Exception:
        if ui:
            ui.messageBox("Error:\n" + traceback.format_exc())
