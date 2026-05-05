#!/usr/bin/env python3
"""
Laser-cut centerline generator for the six Fig. S5 paper motifs.

Source geometry:
Rafsanjani and Pasini, "Bistable auxetic mechanical metamaterials inspired
by ancient geometric motifs", Extreme Mechanics Letters 9 (2016) 291-296,
Supplementary Material Fig. S5.

The circular and parallel motifs are generated from the Fig. S5 building-block
topology, not from only the central rotating unit.  This is important: those
motifs include pinwheel-like cuts that run from the block boundary to the
rotating-unit hinge points.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Iterable, Sequence

from kirigami_geometry_io import (
    EPS,
    Point,
    Polyline,
    reference_template_sheet,
    square_reference_template_sheet,
    write_dxf,
    write_preview_png,
    write_svg,
)


SQRT2 = math.sqrt(2.0)
SQRT3 = math.sqrt(3.0)


def dist(a: Point, b: Point) -> float:
    return (b - a).norm()


def clamp01(x: float) -> float:
    return max(0.0, min(1.0, x))


def theta_square_from_a(a_over_l: float) -> float:
    if not 0.0 <= a_over_l <= 1.0:
        raise ValueError("square tilted a/l must be in [0, 1]")
    return math.acos(a_over_l / SQRT2) - math.pi / 4.0


def theta_triangular_from_a(a_over_l: float) -> float:
    if not 0.0 <= a_over_l <= 1.0:
        raise ValueError("triangular tilted a/l must be in [0, 1]")
    return math.acos(a_over_l / 2.0) - math.pi / 3.0


def square_radius_from_a(a_over_l: float) -> float:
    if not 0.0 <= a_over_l <= (SQRT3 - 1.0) / SQRT2:
        raise ValueError("square circular a/l must be in [0, (sqrt(3)-1)/sqrt(2)]")
    return 0.5 * math.sqrt(1.0 + (1.0 + SQRT2 * a_over_l) ** 2)


def square_a_from_radius(r_over_l: float) -> float:
    if not 1.0 / SQRT2 <= r_over_l <= 1.0:
        raise ValueError("square circular R/l must be in [1/sqrt(2), 1]")
    return (-1.0 + math.sqrt(4.0 * r_over_l * r_over_l - 1.0)) / SQRT2


def triangular_radius_from_a(a_over_l: float) -> float:
    if not 0.0 <= a_over_l <= 1.0:
        raise ValueError("triangular circular a/l must be in [0, 1]")
    return math.sqrt(((2.0 * a_over_l + 1.0) ** 2 + 3.0) / 12.0)


def triangular_a_from_radius(r_over_l: float) -> float:
    if not SQRT3 / 3.0 <= r_over_l <= SQRT3 / 2.0:
        raise ValueError("triangular circular R/l must be in [sqrt(3)/3, sqrt(3)/2]")
    return (-1.0 + math.sqrt(12.0 * r_over_l * r_over_l - 3.0)) / 2.0


def square_width_from_a(a_over_l: float) -> float:
    if not 0.0 <= a_over_l <= 1.0 / SQRT2:
        raise ValueError("square parallel a/l must be in [0, 1/sqrt(2)]")
    return (1.0 - SQRT2 * a_over_l) / 2.0


def square_a_from_width(w_over_l: float) -> float:
    if not 0.0 <= w_over_l <= 0.5:
        raise ValueError("square parallel w/l must be in [0, 1/2]")
    return (1.0 - 2.0 * w_over_l) / SQRT2


def triangular_width_from_a(a_over_l: float) -> float:
    if not 0.0 <= a_over_l <= 1.0:
        raise ValueError("triangular parallel a/l must be in [0, 1]")
    return (1.0 - a_over_l) / (2.0 * SQRT3)


def triangular_a_from_width(w_over_l: float) -> float:
    if not 0.0 <= w_over_l <= SQRT3 / 6.0:
        raise ValueError("triangular parallel w/l must be in [0, sqrt(3)/6]")
    return 1.0 - 2.0 * SQRT3 * w_over_l


def circle_point(center: Point, radius: float, angle: float) -> Point:
    return Point(center.x + radius * math.cos(angle), center.y + radius * math.sin(angle))


def angle(center: Point, point: Point) -> float:
    return math.atan2(point.y - center.y, point.x - center.x)


def shortest_delta(a0: float, a1: float) -> float:
    return (a1 - a0 + math.pi) % (2.0 * math.pi) - math.pi


def arc_polyline(center: Point, radius: float, start: Point, end: Point, segments_per_pi: int) -> Polyline:
    a0 = angle(center, start)
    delta = shortest_delta(a0, angle(center, end))
    n = max(3, int(math.ceil(abs(delta) / (math.pi / max(4, segments_per_pi)))))
    return Polyline(tuple(circle_point(center, radius, a0 + delta * i / n) for i in range(n + 1)))


def arc_chain(center: Point, radius: float, waypoints: Sequence[Point], segments_per_pi: int) -> Polyline:
    points: list[Point] = []
    for start, end in zip(waypoints[:-1], waypoints[1:]):
        segment = arc_polyline(center, radius, start, end, segments_per_pi).points
        if points:
            points.extend(segment[1:])
        else:
            points.extend(segment)
    return Polyline(tuple(points))


def shift_on_circle_by_chord(center: Point, radius: float, point: Point, chord: float, sign: float) -> Point:
    if chord <= EPS:
        return point
    if chord >= 2.0 * radius:
        raise ValueError("circular hinge gap is too large for the selected radius")
    delta = 2.0 * math.asin(chord / (2.0 * radius))
    return circle_point(center, radius, angle(center, point) + sign * delta)


def shift_on_circle_toward(center: Point, radius: float, point: Point, toward: Point, chord: float) -> Point:
    if chord <= EPS:
        return point
    if chord >= 2.0 * radius:
        raise ValueError("circular hinge gap is too large for the selected radius")
    delta = 2.0 * math.asin(chord / (2.0 * radius))
    a0 = angle(center, point)
    candidates = (
        circle_point(center, radius, a0 + delta),
        circle_point(center, radius, a0 - delta),
    )
    return min(candidates, key=lambda p: dist(p, toward))


def shift_segment_toward(point: Point, toward: Point, gap: float) -> Point:
    if gap <= EPS:
        return point
    return point + (toward - point).unit() * gap


def point_at_distance(points: Sequence[Point], s: float) -> Point:
    if s <= 0.0:
        return points[0]
    remaining = s
    for p0, p1 in zip(points[:-1], points[1:]):
        length = dist(p0, p1)
        if remaining <= length + EPS:
            return p0 + (p1 - p0) * clamp01(remaining / length)
        remaining -= length
    return points[-1]


def trim_poly(polyline: Polyline, trim_start: float, trim_end: float) -> Polyline | None:
    total = polyline.length()
    if total <= trim_start + trim_end + EPS:
        return None
    pts = polyline.points
    start = point_at_distance(pts, trim_start)
    end = point_at_distance(pts, total - trim_end)
    kept = [start]
    walked = 0.0
    for p0, p1 in zip(pts[:-1], pts[1:]):
        seg = dist(p0, p1)
        next_walked = walked + seg
        if next_walked > trim_start + EPS and walked < total - trim_end - EPS:
            if walked >= trim_start - EPS and next_walked <= total - trim_end + EPS:
                kept.append(p1)
        walked = next_walked
    if dist(kept[-1], end) > EPS:
        kept.append(end)
    if len(kept) < 2 or dist(kept[0], kept[-1]) <= EPS:
        return None
    return Polyline(tuple(kept))


def trimmed_segment(p0: Point, p1: Point, trim: float) -> Polyline | None:
    return trim_poly(Polyline((p0, p1)), trim, trim)


def add_trimmed(out: list[Polyline], poly: Polyline, trim: float) -> None:
    kept = trim_poly(poly, trim, trim)
    if kept is not None:
        out.append(kept)


def add_trimmed_asym(out: list[Polyline], poly: Polyline, trim_start: float, trim_end: float) -> None:
    kept = trim_poly(poly, trim_start, trim_end)
    if kept is not None:
        out.append(kept)


def square_parallel_block(origin: Point, l: float, w_over_l: float, t: float) -> list[Polyline]:
    w = w_over_l * l
    x0, y0 = origin.x, origin.y
    pb = Point(x0 + 0.5 * l, y0 + w)
    pr = Point(x0 + l - w, y0 + 0.5 * l)
    pt = Point(x0 + 0.5 * l, y0 + l - w)
    pl = Point(x0 + w, y0 + 0.5 * l)
    # At each vertex, the outer branch and one diagonal remain connected.
    # The other diagonal is shifted toward the cell center by the perpendicular
    # clearance t. This makes t the horizontal/vertical gap shown in Fig. S5
    # finite-hinge sketches, not a distance along a slanted cut.
    b_shift = Point(pb.x, pb.y + t)
    r_shift = Point(pr.x - t, pr.y)
    t_shift = Point(pt.x, pt.y - t)
    l_shift = Point(pl.x + t, pl.y)
    return [
        Polyline((Point(pl.x, y0 + l), pl, b_shift)),
        Polyline((Point(x0, pb.y), pb, r_shift)),
        Polyline((Point(pr.x, y0), pr, t_shift)),
        Polyline((Point(x0 + l, pt.y), pt, l_shift)),
    ]


def mirror_square_cell_x(polylines: Sequence[Polyline], x_left: float, l: float) -> list[Polyline]:
    x_right = x_left + l
    return [
        Polyline(tuple(Point(x_right - (p.x - x_left), p.y) for p in poly.points), poly.layer)
        for poly in polylines
    ]


def square_circular_block(origin: Point, l: float, r_over_l: float, t: float, arc_segments: int) -> list[Polyline]:
    r = r_over_l * l
    if r < 0.5 * l:
        raise ValueError("square circular radius is too small")
    h = math.sqrt(max(0.0, r * r - 0.25 * l * l))
    x0, y0 = origin.x, origin.y
    bl = Point(x0, y0)
    br = Point(x0 + l, y0)
    tr = Point(x0 + l, y0 + l)
    tl = Point(x0, y0 + l)

    p_bottom = Point(x0 + 0.5 * l, y0 + l - h)
    p_right = Point(x0 + h, y0 + 0.5 * l)
    p_top = Point(x0 + 0.5 * l, y0 + h)
    p_left = Point(x0 + l - h, y0 + 0.5 * l)

    # Finite hinge treatment matching the square-parallel convention:
    # one arc passes exactly through the ideal rotating-unit corner, while
    # the neighboring arc endpoint is shifted along its circle so the
    # chord gap from that corner is t.
    p_bottom_gap = shift_on_circle_by_chord(tr, r, p_bottom, t, -1.0)
    p_right_gap = shift_on_circle_by_chord(tl, r, p_right, t, -1.0)
    p_top_gap = shift_on_circle_by_chord(bl, r, p_top, t, -1.0)
    p_left_gap = shift_on_circle_by_chord(br, r, p_left, t, -1.0)

    return [
        arc_chain(tl, r, (Point(x0, y0 + l - r), p_bottom, p_right_gap), arc_segments),
        arc_chain(bl, r, (Point(x0 + r, y0), p_right, p_top_gap), arc_segments),
        arc_chain(br, r, (p_left_gap, p_top, Point(x0 + l, y0 + r)), arc_segments),
        arc_chain(tr, r, (p_bottom_gap, p_left, Point(x0 + l - r, y0 + l)), arc_segments),
    ]


def triangle_frame(origin: Point, l: float, col: int, row: int, up: bool) -> tuple[Point, Point, Point]:
    e1 = Point(l, 0.0)
    e2 = Point(0.5 * l, 0.5 * SQRT3 * l)
    o = origin + e1 * col + e2 * row
    if up:
        return o, o + e1, o + e2
    return o + e1, o + e1 + e2, o + e2


def tri_local(A: Point, B: Point, x: float, y: float) -> Point:
    u = (B - A).unit()
    v = u.left_normal()
    return A + u * x + v * y


def line_intersection(p1: Point, d1: Point, p2: Point, d2: Point) -> Point:
    den = d1.cross(d2)
    if abs(den) < EPS:
        raise ValueError("parallel lines do not intersect")
    return p1 + d1 * ((p2 - p1).cross(d2) / den)


def triangular_parallel_block(triangle: tuple[Point, Point, Point], l: float, w_over_l: float, t: float) -> list[Polyline]:
    A, B, C = triangle
    w = w_over_l * l
    h = 0.5 * SQRT3 * l
    p_left = tri_local(A, B, SQRT3 * w, w)
    p_right = tri_local(A, B, l - SQRT3 * w, w)
    p_top = tri_local(A, B, 0.5 * l, h - 2.0 * w)
    e_left = tri_local(A, B, w / SQRT3, w)
    e_bottom = tri_local(A, B, l - 2.0 * w / SQRT3, 0.0)

    side_dir = (C - B).unit()
    e_side = line_intersection(p_left, p_top - p_left, B, side_dir)

    p_right_gap = shift_segment_toward(p_right, p_left, t)
    p_top_gap = shift_segment_toward(p_top, p_right, t)
    p_left_gap = shift_segment_toward(p_left, p_top, t)
    return [
        Polyline((e_left, p_left, p_right_gap)),
        Polyline((e_bottom, p_right, p_top_gap)),
        Polyline((e_side, p_top, p_left_gap)),
    ]


def triangular_circular_block(
    triangle: tuple[Point, Point, Point],
    l: float,
    r_over_l: float,
    t: float,
    arc_segments: int,
) -> list[Polyline]:
    A, B, C = triangle
    r = r_over_l * l
    a_over_l = (-1.0 + math.sqrt(max(0.0, 12.0 * r_over_l * r_over_l - 3.0))) / 2.0
    a = a_over_l * l
    h = 0.5 * SQRT3 * l
    half_chord = 0.5 * a
    drop = math.sqrt(max(0.0, r * r - half_chord * half_chord))
    y_base = h - drop

    p_left = tri_local(A, B, 0.5 * l - half_chord, y_base)
    p_right = tri_local(A, B, 0.5 * l + half_chord, y_base)
    p_top = tri_local(A, B, 0.5 * l, y_base + 0.5 * SQRT3 * a)
    e_left = A + (C - A) * (1.0 - r_over_l)
    e_bottom = tri_local(A, B, r, 0.0)
    e_side = B + (C - B) * r_over_l

    # Same topology as triangular parallel: three connected arc chains.  Each
    # chain passes exactly through one rotating-unit corner and stops t away
    # from the next corner along the same circle.
    p_right_gap = shift_on_circle_toward(C, r, p_right, p_left, t)
    p_top_gap = shift_on_circle_toward(A, r, p_top, p_right, t)
    p_left_gap = shift_on_circle_toward(B, r, p_left, p_top, t)
    return [
        arc_chain(C, r, (e_left, p_left, p_right_gap), arc_segments),
        arc_chain(A, r, (e_bottom, p_right, p_top_gap), arc_segments),
        arc_chain(B, r, (e_side, p_top, p_left_gap), arc_segments),
    ]


def clip_segment_to_rect(p0: Point, p1: Point, width: float, height: float) -> tuple[Point, Point] | None:
    dx = p1.x - p0.x
    dy = p1.y - p0.y
    u0, u1 = 0.0, 1.0
    for p, q in ((-dx, p0.x), (dx, width - p0.x), (-dy, p0.y), (dy, height - p0.y)):
        if abs(p) < EPS:
            if q < -EPS:
                return None
            continue
        r = q / p
        if p < 0.0:
            u0 = max(u0, r)
        else:
            u1 = min(u1, r)
        if u0 > u1:
            return None
    a = Point(p0.x + u0 * dx, p0.y + u0 * dy)
    b = Point(p0.x + u1 * dx, p0.y + u1 * dy)
    if dist(a, b) <= EPS:
        return None
    return a, b


def clip_polylines_to_rect(polylines: Sequence[Polyline], width: float, height: float) -> list[Polyline]:
    out: list[Polyline] = []
    for poly in polylines:
        if all(-EPS <= p.x <= width + EPS and -EPS <= p.y <= height + EPS for p in poly.points):
            out.append(poly)
            continue
        for p0, p1 in zip(poly.points[:-1], poly.points[1:]):
            clipped = clip_segment_to_rect(p0, p1, width, height)
            if clipped is not None:
                out.append(Polyline(clipped, poly.layer))
    return out


def resolve_handle_dimensions(
    width: float,
    height: float,
    handle_width: float | None,
    handle_height: float | None,
) -> tuple[float, float]:
    short_side = min(width, height)
    resolved_width = 0.25 * short_side if handle_width is None else handle_width
    resolved_height = 0.15 * short_side if handle_height is None else handle_height
    if resolved_width <= 0.0:
        raise ValueError("--handle-width must be positive")
    if resolved_height <= 0.0:
        raise ValueError("--handle-height must be positive")
    return resolved_width, resolved_height


def selected_handle_sides(side: str) -> set[str]:
    mapping = {
        "top": {"top"},
        "bottom": {"bottom"},
        "left": {"left"},
        "right": {"right"},
        "both": {"top", "bottom"},
        "top-bottom": {"top", "bottom"},
        "left-right": {"left", "right"},
        "all": {"top", "bottom", "left", "right"},
    }
    if side not in mapping:
        raise ValueError("--handle-side must be one of: top, bottom, left, right, top-bottom, left-right, both, all")
    return mapping[side]


def outline_with_optional_handle(args: argparse.Namespace) -> tuple[list[Polyline], dict[str, float | str]]:
    width = args.width
    height = args.height
    if not getattr(args, "handle", False):
        return [
            Polyline(
                (Point(0.0, 0.0), Point(width, 0.0), Point(width, height), Point(0.0, height), Point(0.0, 0.0)),
                "OUTLINE",
            )
        ], {}

    handle_width, handle_height = resolve_handle_dimensions(
        width,
        height,
        getattr(args, "handle_width", None),
        getattr(args, "handle_height", None),
    )
    side = getattr(args, "handle_side", "both")
    sides = selected_handle_sides(side)
    if sides & {"top", "bottom"} and handle_width >= width:
        raise ValueError("--handle-width must be smaller than the sheet width for top/bottom handles")
    if sides & {"left", "right"} and handle_width >= height:
        raise ValueError("--handle-width must be smaller than the sheet height for left/right handles")

    cx = 0.5 * width
    x_left = cx - 0.5 * handle_width
    x_right = cx + 0.5 * handle_width
    cy = 0.5 * height
    y_bottom = cy - 0.5 * handle_width
    y_top = cy + 0.5 * handle_width

    points: list[Point] = []

    def add(point: Point) -> None:
        if not points or dist(points[-1], point) > EPS:
            points.append(point)

    add(Point(0.0, 0.0))
    add(Point(x_left, 0.0))
    if "bottom" in sides:
        add(Point(x_left, -handle_height))
        add(Point(x_right, -handle_height))
        add(Point(x_right, 0.0))
    add(Point(width, 0.0))
    add(Point(width, y_bottom))
    if "right" in sides:
        add(Point(width + handle_height, y_bottom))
        add(Point(width + handle_height, y_top))
        add(Point(width, y_top))
    add(Point(width, height))
    add(Point(x_right, height))
    if "top" in sides:
        add(Point(x_right, height + handle_height))
        add(Point(x_left, height + handle_height))
        add(Point(x_left, height))
    add(Point(0.0, height))
    add(Point(0.0, y_top))
    if "left" in sides:
        add(Point(-handle_height, y_top))
        add(Point(-handle_height, y_bottom))
        add(Point(0.0, y_bottom))
    add(Point(0.0, 0.0))
    return [Polyline(tuple(points), "OUTLINE")], {
        "handle": side,
        "handle_width": handle_width,
        "handle_height": handle_height,
    }


def triangular_checkerboard_sheet(
    motif_paths: Sequence[Polyline],
    width: float,
    height: float,
    l: float,
    offset_x: float,
    offset_y: float,
) -> list[Polyline]:
    Matrix = tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]

    def matmul(a: Matrix, b: Matrix) -> Matrix:
        return tuple(
            tuple(sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3))
            for i in range(3)
        )  # type: ignore[return-value]

    def translate(dx: float, dy: float) -> Matrix:
        return ((1.0, 0.0, dx), (0.0, 1.0, dy), (0.0, 0.0, 1.0))

    def scale(sx: float, sy: float) -> Matrix:
        return ((sx, 0.0, 0.0), (0.0, sy, 0.0), (0.0, 0.0, 1.0))

    def rotate(angle_deg: float, about_x: float = 0.0, about_y: float = 0.0) -> Matrix:
        a = math.radians(angle_deg)
        c = math.cos(a)
        s = math.sin(a)
        r = ((c, -s, 0.0), (s, c, 0.0), (0.0, 0.0, 1.0))
        return matmul(translate(about_x, about_y), matmul(r, translate(-about_x, -about_y)))

    def combine(transforms: Sequence[Matrix]) -> Matrix:
        matrix: Matrix = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
        for transform in transforms:
            matrix = matmul(matrix, transform)
        return matrix

    def transform_point(matrix: Matrix, p: Point) -> Point:
        return Point(
            matrix[0][0] * p.x + matrix[0][1] * p.y + matrix[0][2],
            matrix[1][0] * p.x + matrix[1][1] * p.y + matrix[1][2],
        )

    def to_cad(p: Point) -> Point:
        return Point(p.x, height - p.y)

    cx = 0.5 * l
    cy = SQRT3 * l / 6.0
    x_step = 0.5 * l
    y_step = 0.5 * SQRT3 * l
    cols = max(1, int(math.ceil(width / x_step)))
    rows = max(1, int(math.ceil(height / y_step)))
    outer = translate(offset_x - 0.25 * l, offset_y)

    paths: list[Polyline] = []
    for col in range(cols):
        for row in range(rows):
            cell_transforms: list[Matrix] = [outer, translate(col * x_step, row * y_step)]
            if (col + row) % 2 == 1:
                cell_transforms.extend(
                    [
                        translate(l, cy),
                        scale(-1.0, 1.0),
                        rotate(180.0, cx, cy),
                    ]
                )
            matrix = combine(cell_transforms)
            transformed = [
                Polyline(tuple(transform_point(matrix, p) for p in poly.points), poly.layer)
                for poly in motif_paths
            ]
            for clipped in clip_polylines_to_rect(transformed, width, height):
                paths.append(Polyline(tuple(to_cad(p) for p in clipped.points), clipped.layer))
    return paths


def generate_sheet(args: argparse.Namespace) -> tuple[list[Polyline], list[Polyline], dict[str, float | str]]:
    width = args.width
    height = args.height
    l = args.cell_size
    a_over_l = args.a_over_l
    t = args.t
    family = args.family
    motif = args.motif
    info: dict[str, float | str] = {"family": family, "motif": motif, "l": l, "t": t, "a/l": a_over_l}

    if motif == "tilted":
        if args.theta_deg is not None:
            theta = math.radians(args.theta_deg)
            if family == "square":
                if not 0.0 <= theta <= math.pi / 4.0:
                    raise ValueError("square tilted theta must be in [0, 45] degrees")
                info["a/l"] = math.cos(theta) - math.sin(theta)
            else:
                if not 0.0 <= theta <= math.pi / 6.0:
                    raise ValueError("triangular tilted theta must be in [0, 30] degrees")
                info["a/l"] = math.cos(theta) - SQRT3 * math.sin(theta)
        elif family == "square":
            theta = theta_square_from_a(a_over_l)
        else:
            theta = theta_triangular_from_a(a_over_l)
        info["theta_deg"] = math.degrees(theta)
        if family == "square":
            return (*square_reference_template_sheet(width, height, l, theta, t), info)
        return (*reference_template_sheet(width, height, l, theta, t), info)

    paths: list[Polyline] = []
    if family == "square":
        if motif == "circular":
            if args.radius_over_l is not None:
                r_over_l = args.radius_over_l
                info["a/l"] = square_a_from_radius(r_over_l)
            else:
                r_over_l = square_radius_from_a(a_over_l)
            info["R/l"] = r_over_l
        else:
            if args.width_over_l is not None:
                w_over_l = args.width_over_l
                info["a/l"] = square_a_from_width(w_over_l)
            else:
                w_over_l = square_width_from_a(a_over_l)
            info["w/l"] = w_over_l
        col_min = math.floor(-args.offset_x / l) - 2
        col_max = math.ceil((width - args.offset_x) / l) + 2
        row_min = math.floor(-args.offset_y / l) - 2
        row_max = math.ceil((height - args.offset_y) / l) + 2
        for row in range(row_min, row_max + 1):
            for col in range(col_min, col_max + 1):
                o = Point(args.offset_x + col * l, args.offset_y + row * l)
                if motif == "circular":
                    block_paths = square_circular_block(o, l, float(info["R/l"]), t, args.arc_segments)
                    if (col + row) % 2 == 1:
                        block_paths = mirror_square_cell_x(block_paths, o.x, l)
                else:
                    block_paths = square_parallel_block(o, l, float(info["w/l"]), t)
                    if (col + row) % 2 == 1:
                        block_paths = mirror_square_cell_x(block_paths, o.x, l)
                paths.extend(block_paths)
        if motif in {"circular", "parallel"}:
            info["cell_orientation"] = "checkerboard: normal when (col+row) even, mirrored-x when odd"
    else:
        if motif == "circular":
            if args.radius_over_l is not None:
                r_over_l = args.radius_over_l
                info["a/l"] = triangular_a_from_radius(r_over_l)
            else:
                r_over_l = triangular_radius_from_a(a_over_l)
            info["R/l"] = r_over_l
        else:
            if args.width_over_l is not None:
                w_over_l = args.width_over_l
                info["a/l"] = triangular_a_from_width(w_over_l)
            else:
                w_over_l = triangular_width_from_a(a_over_l)
            info["w/l"] = w_over_l

        canonical_triangle = (Point(0.0, 0.0), Point(l, 0.0), Point(0.5 * l, 0.5 * SQRT3 * l))
        if motif == "circular":
            motif_paths = triangular_circular_block(
                canonical_triangle,
                l,
                float(info["R/l"]),
                t,
                args.arc_segments,
            )
        else:
            motif_paths = triangular_parallel_block(canonical_triangle, l, float(info["w/l"]), t)
        paths.extend(triangular_checkerboard_sheet(motif_paths, width, height, l, args.offset_x, args.offset_y))
        info["cell_orientation"] = "tilted-triangle checkerboard placement: normal when (col+row) even, mirrored-x/rotated when odd"

    outline, outline_info = outline_with_optional_handle(args)
    info.update(outline_info)
    return clip_polylines_to_rect(paths, width, height), outline, info


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--family", choices=("square", "triangular"), required=True)
    parser.add_argument("--motif", choices=("tilted", "circular", "parallel"), required=True)
    parser.add_argument("--width", type=float, required=True, help="sheet width in mm")
    parser.add_argument("--height", type=float, required=True, help="sheet height in mm")
    parser.add_argument("--cell-size", type=float, required=True, help="paper building-block side length l in mm")
    parser.add_argument("--t", type=float, default=1.0, help="finite cut-tip trim / ligament scale in mm")
    parser.add_argument("--a-over-l", type=float, default=0.5, help="rotating-unit length ratio a/l")
    parser.add_argument("--theta-deg", type=float, default=None, help="tilted motif angle override")
    parser.add_argument("--radius-over-l", type=float, default=None, help="circular motif R/l override")
    parser.add_argument("--width-over-l", type=float, default=None, help="parallel motif w/l override")
    parser.add_argument("--offset-x", type=float, default=0.0)
    parser.add_argument("--offset-y", type=float, default=0.0)
    parser.add_argument("--arc-segments", type=int, default=24)
    parser.add_argument("--cut-width", type=float, default=0.2)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output"),
        help="folder used for automatic output names when --svg/--dxf/--preview-png are omitted",
    )
    parser.add_argument("--svg", type=Path, default=None)
    parser.add_argument("--dxf", type=Path, default=None)
    parser.add_argument("--preview-png", type=Path, default=None)
    parser.add_argument("--preview-scale", type=float, default=8.0)
    parser.add_argument("--preview-line-width", type=int, default=3)
    parser.add_argument("--handle", action="store_true", help="add centered rectangular handle tab(s) to the contour")
    parser.add_argument(
        "--handle-side",
        choices=("top", "bottom", "left", "right", "both", "top-bottom", "left-right", "all"),
        default="both",
        help="which side receives a handle tab when --handle is used",
    )
    parser.add_argument("--handle-width", type=float, default=None, help="handle tab width in mm; default is 25%% of the short sheet side")
    parser.add_argument("--handle-height", type=float, default=None, help="handle tab protrusion in mm; default is 15%% of the short sheet side")
    parser.add_argument("--no-outline", action="store_true")
    return parser.parse_args()


def filename_token(value: float | str) -> str:
    if isinstance(value, float):
        return (
            f"{value:.6g}"
            .replace("-", "m")
            .replace("+", "")
            .replace(".", "p")
            .replace("/", "over")
            .replace(" ", "")
        )
    else:
        text = str(value)
    return (
        text.replace("-", "")
        .replace("+", "")
        .replace(".", "p")
        .replace("/", "over")
        .replace(" ", "")
    )


def automatic_output_stem(args: argparse.Namespace, info: dict[str, float | str]) -> str:
    parts = [
        args.family,
        args.motif,
        f"W{filename_token(args.width)}",
        f"H{filename_token(args.height)}",
        f"l{filename_token(args.cell_size)}",
        f"t{filename_token(args.t)}",
        f"aol{filename_token(info.get('a/l', args.a_over_l))}",
    ]
    if "theta_deg" in info:
        parts.append(f"theta{filename_token(info['theta_deg'])}deg")
    if "R/l" in info:
        parts.append(f"Rol{filename_token(info['R/l'])}")
    if "w/l" in info:
        parts.append(f"wol{filename_token(info['w/l'])}")
    if abs(args.offset_x) > EPS or abs(args.offset_y) > EPS:
        parts.extend([f"ox{filename_token(args.offset_x)}", f"oy{filename_token(args.offset_y)}"])
    if getattr(args, "handle", False):
        parts.append(f"handle{filename_token(info.get('handle', args.handle_side))}")
        parts.append(f"hw{filename_token(info.get('handle_width', 0.0))}")
        parts.append(f"hh{filename_token(info.get('handle_height', 0.0))}")
    if args.no_outline:
        parts.append("nooutline")
    return "_".join(parts)


def fill_automatic_outputs(args: argparse.Namespace, info: dict[str, float | str]) -> None:
    if args.svg is not None or args.dxf is not None or args.preview_png is not None:
        return
    stem = automatic_output_stem(args, info)
    args.svg = args.output_dir / f"{stem}.svg"
    args.dxf = args.output_dir / f"{stem}.dxf"
    args.preview_png = args.output_dir / f"{stem}.png"


def main() -> None:
    args = parse_args()
    cut_paths, outline_paths, info = generate_sheet(args)
    if args.handle:
        outline_paths, outline_info = outline_with_optional_handle(args)
        info.update(outline_info)
    if args.no_outline:
        outline_paths = []
    fill_automatic_outputs(args, info)
    for output in (args.svg, args.dxf, args.preview_png):
        if output is not None:
            output.parent.mkdir(parents=True, exist_ok=True)
    if args.svg is not None:
        write_svg(args.svg, cut_paths, outline_paths, args.cut_width, margin=0.0)
    if args.dxf is not None:
        write_dxf(args.dxf, cut_paths, outline_paths)
    if args.preview_png is not None:
        write_preview_png(args.preview_png, cut_paths, outline_paths, scale=args.preview_scale, line_width=args.preview_line_width)

    print("Fig. S5 paper-motif laser pattern")
    for key, value in info.items():
        if isinstance(value, float):
            print(f"  {key}: {value:.8g}")
        else:
            print(f"  {key}: {value}")
    print(f"  sheet: {args.width:.8g} x {args.height:.8g} mm")
    print(f"  cut polylines: {len(cut_paths)}")
    if args.svg is not None:
        print(f"  SVG: {args.svg}")
    if args.dxf is not None:
        print(f"  DXF: {args.dxf}")
    if args.preview_png is not None:
        print(f"  preview: {args.preview_png}")


if __name__ == "__main__":
    main()
