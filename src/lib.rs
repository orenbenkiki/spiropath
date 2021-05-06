// Copyright (C) 2021 Oren Ben-Kiki
//
// This program is free software: you can redistribute it and/or modify it under the terms of the
// GNU Affero General Public License as published by the Free Software Foundation, either version 3
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
// even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License along with this program.
// If not, see <https://www.gnu.org/licenses/>.

//! Generate Spirograph-like paths.

#[macro_use]
extern crate assert_float_eq;

use ordered_float::OrderedFloat;
use std::cmp::Ordering;
use std::f64::consts::PI;
use std::fs::read_to_string;
use std::io::Write;
#[cfg(test)]
use std::str::from_utf8;
use svg2polylines::parse as parse_svg;
use svg2polylines::CoordinatePair as Point;
use svg2polylines::Polyline;

fn sign(ordering: Ordering) -> isize {
    match ordering {
        Ordering::Less => -1,
        Ordering::Equal => 0,
        Ordering::Greater => 1,
    }
}

/// Test whether a polyline has a clockwise orientation.
pub fn is_clockwise(polyline: &[Point]) -> bool {
    let mut sum: isize = 0;

    let left_index = (0..polyline.len())
        .min_by_key(|index| OrderedFloat(polyline[*index].x))
        .unwrap();
    let left_y_before =
        OrderedFloat(polyline[(left_index + polyline.len() - 1) % polyline.len()].y);
    let left_y_after = OrderedFloat(polyline[(left_index + 1) % polyline.len()].y);
    sum -= sign(left_y_before.cmp(&left_y_after));

    let right_index = (0..polyline.len())
        .max_by_key(|index| OrderedFloat(polyline[*index].x))
        .unwrap();
    let right_y_before =
        OrderedFloat(polyline[(right_index + polyline.len() - 1) % polyline.len()].y);
    let right_y_after = OrderedFloat(polyline[(right_index + 1) % polyline.len()].y);
    sum += sign(right_y_before.cmp(&right_y_after));

    let top_index = (0..polyline.len())
        .max_by_key(|index| OrderedFloat(polyline[*index].y))
        .unwrap();
    let top_x_before = OrderedFloat(polyline[(top_index + polyline.len() - 1) % polyline.len()].x);
    let top_x_after = OrderedFloat(polyline[(top_index + 1) % polyline.len()].x);
    sum -= sign(top_x_before.cmp(&top_x_after));

    let bottom_index = (0..polyline.len())
        .min_by_key(|index| OrderedFloat(polyline[*index].y))
        .unwrap();
    let bottom_x_before =
        OrderedFloat(polyline[(bottom_index + polyline.len() - 1) % polyline.len()].x);
    let bottom_x_after = OrderedFloat(polyline[(bottom_index + 1) % polyline.len()].x);
    sum += sign(bottom_x_before.cmp(&bottom_x_after));

    sum >= 0
}

#[cfg(test)]
#[test]
fn test_clockwise() {
    let mut triangle = vec![
        Point { x: 0.0, y: 0.0 },
        Point { x: 1.0, y: 0.0 },
        Point { x: 1.0, y: 1.0 },
    ];
    assert!(!is_clockwise(&triangle));
    triangle.reverse();
    assert!(is_clockwise(&triangle));
}

/// Return the distance between two points.
pub fn distance(left_point: Point, right_point: Point) -> f64 {
    let delta_x = left_point.x - right_point.x;
    let delta_y = left_point.y - right_point.y;
    delta_x.hypot(delta_y)
}

/// Measure the lengths of each of the segments of the polyline.
pub fn lengths(polyline: &[Point]) -> Vec<f64> {
    (0..polyline.len())
        .map(|prev_index| {
            let next_index = (prev_index + 1) % polyline.len();
            distance(polyline[prev_index], polyline[next_index])
        })
        .collect()
}

/// Scale a polyline by a factor.
pub fn scale_lengths(lengths: &mut Vec<f64>, factor: f64) {
    for length in lengths.iter_mut() {
        *length *= factor;
    }
}

/// Scale a point by a factor.
pub fn scale_point(point: &mut Point, x_factor: f64, y_factor: f64) {
    point.x *= x_factor;
    point.y *= y_factor;
}

/// Scale a polyline by a factor.
pub fn scale_polyline(polyline: &mut Polyline, x_factor: f64, y_factor: f64) {
    for mut point in polyline.iter_mut() {
        scale_point(&mut point, x_factor, y_factor);
    }
}

#[cfg(test)]
#[test]
fn test_scale() {
    let mut square = vec![
        Point { x: 0.0, y: 0.0 },
        Point { x: 1.0, y: 0.0 },
        Point { x: 1.0, y: 1.0 },
        Point { x: 0.0, y: 1.0 },
    ];
    assert!(lengths(&square).iter().sum::<f64>() == 4.0);
    scale_polyline(&mut square, 2.0, 3.0);
    assert!(lengths(&square).iter().sum::<f64>() == 10.0);
}

/// Return the direction (length 1) vector between two points (as a point).
pub fn direction(from_point: Point, to_point: Point) -> Point {
    let delta_x = to_point.x - from_point.x;
    let delta_y = to_point.y - from_point.y;
    let distance = delta_x.hypot(delta_y);
    Point {
        x: delta_x / distance,
        y: delta_y / distance,
    }
}

/// Return the direction of the path at each of its points.
pub fn directions(polyline: &[Point]) -> Polyline {
    let mut result = vec![Point { x: 0.0, y: 0.0 }; polyline.len()];
    for from_index in 0..polyline.len() {
        let to_index = (from_index + 1) % polyline.len();
        result[from_index] = direction(polyline[from_index], polyline[to_index]);
    }
    result
}

#[cfg(test)]
fn assert_point(point: Point, x: f64, y: f64) {
    assert_float_absolute_eq!(point.x, x, 1e-6);
    assert_float_absolute_eq!(point.y, y, 1e-6);
}

#[cfg(test)]
#[test]
fn test_directions() {
    let triangle = vec![
        Point { x: 0.0, y: 0.0 },
        Point { x: 1.0, y: 0.0 },
        Point { x: 1.0, y: 1.0 },
    ];
    let triangle_directions = directions(&triangle);

    assert_point(triangle_directions[0], 1.0, 0.0);
    assert_point(triangle_directions[1], 0.0, 1.0);
    assert_point(triangle_directions[2], -0.5_f64.sqrt(), -0.5_f64.sqrt());
}

/// Find the index of the top-most point of a polyline (after rotating by an angle).
pub fn find_top_most_point_index(polyline: &[Point], degrees: f64) -> usize {
    let radians = degrees * PI / 180.0;
    let sine = radians.sin();
    let cosine = radians.cos();
    (0..polyline.len())
        .max_by_key(|index| OrderedFloat(polyline[*index].y * cosine + polyline[*index].x * sine))
        .unwrap()
}

#[cfg(test)]
#[test]
fn test_find_top_most_point_index() {
    let square = vec![
        Point { x: 0.0, y: 0.0 },
        Point { x: 0.0, y: 1.0 },
        Point { x: 1.0, y: 1.0 },
        Point { x: 1.0, y: 0.0 },
    ];
    assert!(find_top_most_point_index(&square, 45.0) == 2);
    assert!(find_top_most_point_index(&square, 135.0) == 3);
    assert!(find_top_most_point_index(&square, -45.0) == 1);
}

// BEGIN NOT TESTED

pub fn prune_repeated_points(polyline: &mut Polyline, tolerance: f64) {
    let mut kept_count = 1;
    for next_index in 1..polyline.len() {
        if distance(polyline[kept_count - 1], polyline[next_index]) <= tolerance {
            continue;
        }
        if kept_count < next_index {
            polyline[kept_count] = polyline[next_index];
            kept_count += 1;
        }
    }

    if kept_count > 1 && distance(polyline[0], polyline[kept_count - 1]) <= tolerance {
        kept_count -= 1;
    }

    polyline.truncate(kept_count);
}

/// Load a path from a file and convert it to a polyline.
pub fn load_polyline_from_svg_file(path: &str, tolerance: f64) -> Polyline {
    let string = read_to_string(path).unwrap_or_else(|_| panic!("reading file: {}", path));

    let mut polylines = parse_svg(&string).unwrap_or_else(|_| panic!("parsing file: {}", path));

    if polylines.is_empty() {
        panic!("no SVG paths in file: {}", path);
    }

    if polylines.len() > 1 {
        panic!("too many SVG paths in file: {}", path);
    }

    let mut polyline = polylines.pop().unwrap();

    if polyline.len() < 2 {
        panic!("too few points in SVG path in file: {}", path);
    }

    prune_repeated_points(&mut polyline, tolerance);

    if polyline.len() == 1 {
        panic!("zero length path in file: {}", path);
    }

    polyline
}

// END NOT TESTED

fn print_svg_polyline(polyline: &[Point], output: &mut dyn Write) {
    writeln!(output, "<path d='").unwrap();
    let mut command = "M";
    for point in polyline {
        writeln!(output, "{} {}pt {}pt", command, point.x, point.y).unwrap();
        command = "L";
    }
    writeln!(output, "Z").unwrap();
    writeln!(output, "'/>").unwrap();
}

/// Print a vector of paths as an SVG file.
pub fn print_svg_polylines(polylines: &[Polyline], output: &mut dyn Write) {
    let maximal_point = maximal_coordinates(polylines);
    writeln!(
        output,
        "<svg width='{}pt', height='{}pt', xmlns='http://www.w3.org/2000/svg'>",
        maximal_point.x, maximal_point.y
    )
    .unwrap();
    for polyline in polylines {
        print_svg_polyline(polyline, output);
    }
    writeln!(output, "</svg>").unwrap();
}

/// Location of the rotating shape relative to the stationary shape.
#[derive(Clone, Copy, Debug)]
pub enum Location {
    /// The rotating shape is outside the stationary shape.
    Outside,

    /// The rotating shape in inside the stationary shape.
    Inside,
}

/// How to scale an axis of the output.
#[derive(Clone, Copy, Debug)]
pub enum Scale {
    /// Scale to ensure a given total size.
    Size(f64),

    /// Scale by a fixed factor.
    Factor(f64),

    /// Use the same scale factor as the other axis.
    Same,
}

/// Options for controlling a spiropath.
#[derive(Debug)]
pub struct Options {
    /// Location of the rotating shape relative to the stationary shape.
    pub location: Location,

    /// The number of teeth on the stationary shape.
    pub stationary_teeth: usize,

    /// The number of teeth on the rotation shape.
    pub rotating_teeth: usize,

    /// The tolerance as a fraction of a tooth size.
    pub tolerance: f64,

    /// The point to trace (in the original coordinates of the rotated shape).
    pub traced_point: Point,

    /// The angle of the rotating shape at the start of the stationary shape (in degrees).
    pub initial_rotating_angle: f64,

    /// Whether to include the stationary shape in the output.
    pub include_stationary: bool,

    /// The initial offsets of the rotating shape.
    pub initial_offsets: Vec<f64>,

    /// How to scale the X axis of the output.
    pub x_scale: Scale,

    /// How to scale the Y axis of the output.
    pub y_scale: Scale,
}

impl Options {
    /// Ensure the options are valid.
    pub fn validate(&mut self) {
        if self.stationary_teeth == 0 {
            panic!("stationary-teeth is zero"); // NOT TESTED
        }

        if self.rotating_teeth == 0 {
            panic!("rotating-teeth is zero"); // NOT TESTED
        }

        if self.tolerance <= 0.0 {
            panic!("tolerance {} is not positive", self.tolerance); // NOT TESTED
        }

        if self.initial_offsets.is_empty() {
            self.initial_offsets.push(0.0);
        }
    }
}

/// A point along a polyline where it is in contact with another.
#[derive(Clone, Copy, Debug)]
pub struct ContactPoint {
    /// The index of the point in the polyline.
    pub index: usize,

    /// How far to move along the line segment starting at that point.
    pub offset: f64,
}

/// Find a point along a polyline by its offset (up to some tolerance).
fn find_offset_contact_point(lengths: &[f64], mut offset: f64, mut tolerance: f64) -> ContactPoint {
    tolerance /= 2.0;
    let mut index = 0;

    while offset > tolerance && lengths[index] < offset + tolerance {
        offset -= lengths[index];
        index = (index + 1) % lengths.len();
    }

    if offset <= tolerance {
        offset = 0.0;
    }

    ContactPoint { index, offset }
}

#[cfg(test)]
#[test]
fn test_find_offset_contact_point() {
    let triangle_lengths = vec![1.0, 1.001, 1.0];

    assert!(find_offset_contact_point(&triangle_lengths, 0.0, 0.01).index == 0);
    assert_float_absolute_eq!(
        find_offset_contact_point(&triangle_lengths, 0.0, 0.0).offset,
        0.0,
        1e-6
    );

    assert!(find_offset_contact_point(&triangle_lengths, 1.0, 0.01).index == 1);
    assert_float_absolute_eq!(
        find_offset_contact_point(&triangle_lengths, 1.0, 0.01).offset,
        0.0,
        1e-6
    );

    assert!(find_offset_contact_point(&triangle_lengths, 1.5, 0.01).index == 1);
    assert_float_absolute_eq!(
        find_offset_contact_point(&triangle_lengths, 1.5, 0.01).offset,
        0.5,
        1e-6
    );

    assert!(find_offset_contact_point(&triangle_lengths, 2.0, 0.01).index == 2);
    assert_float_absolute_eq!(
        find_offset_contact_point(&triangle_lengths, 2.0, 0.01).offset,
        0.0,
        1e-6
    );
}

/// Return the minimal coordinates in some paths.
pub fn minimal_coordinates(polylines: &[Polyline]) -> Point {
    let minimal_x = polylines
        .iter()
        .map(|polyline| {
            polyline
                .iter()
                .map(|point| OrderedFloat(point.x))
                .min()
                .unwrap()
        })
        .min()
        .unwrap();

    let minimal_y = polylines
        .iter()
        .map(|polyline| {
            polyline
                .iter()
                .map(|point| OrderedFloat(point.y))
                .min()
                .unwrap()
        })
        .min()
        .unwrap();

    Point {
        x: *minimal_x,
        y: *minimal_y,
    }
}

/// Return the maximal coordinates in some paths.
pub fn maximal_coordinates(polylines: &[Polyline]) -> Point {
    let maximal_x = polylines
        .iter()
        .map(|polyline| {
            polyline
                .iter()
                .map(|point| OrderedFloat(point.x))
                .max()
                .unwrap()
        })
        .max()
        .unwrap();

    let maximal_y = polylines
        .iter()
        .map(|polyline| {
            polyline
                .iter()
                .map(|point| OrderedFloat(point.y))
                .max()
                .unwrap()
        })
        .max()
        .unwrap();

    Point {
        x: *maximal_x,
        y: *maximal_y,
    }
}

/// Scale and move the points to fit in a (0,0) -> (w,h) bounding box.
pub fn transform(polylines: &mut [Polyline], x_scale: Scale, y_scale: Scale) {
    let minimal_point = minimal_coordinates(polylines);
    let maximal_point = maximal_coordinates(polylines);
    let size = Point {
        x: maximal_point.x - minimal_point.x,
        y: maximal_point.y - minimal_point.y,
    };

    let maybe_x_factor = match x_scale {
        Scale::Size(x_size) => Some(x_size / size.x),
        Scale::Factor(factor) => Some(factor),
        Scale::Same => None,
    };

    let maybe_y_factor = match y_scale {
        Scale::Size(y_size) => Some(y_size / size.y),
        Scale::Factor(factor) => Some(factor),
        Scale::Same => None,
    };

    let (x_factor, y_factor) = match (maybe_x_factor, maybe_y_factor) {
        (None, None) => (1.0, 1.0),
        (None, Some(y_factor)) => (y_factor, y_factor),
        (Some(x_factor), None) => (x_factor, x_factor),
        (Some(x_factor), Some(y_factor)) => (x_factor, y_factor),
    };

    for polyline in polylines.iter_mut() {
        for point in polyline.iter_mut() {
            point.x = (point.x - minimal_point.x) * x_factor;
            point.y = (point.y - minimal_point.y) * y_factor;
        }
    }
}

#[cfg(test)]
#[test]
fn test_transform() {
    let mut polylines = vec![vec![
        Point { x: -1.0, y: 0.0 },
        Point { x: 0.0, y: -1.0 },
        Point { x: 1.0, y: 0.0 },
        Point { x: 0.0, y: 1.0 },
    ]];

    transform(&mut polylines, Scale::Same, Scale::Same);
    assert_point(polylines[0][0], 0.0, 1.0);
    assert_point(polylines[0][1], 1.0, 0.0);
    assert_point(polylines[0][2], 2.0, 1.0);
    assert_point(polylines[0][3], 1.0, 2.0);

    transform(&mut polylines, Scale::Factor(2.0), Scale::Same);
    assert_point(polylines[0][0], 0.0, 2.0);
    assert_point(polylines[0][1], 2.0, 0.0);
    assert_point(polylines[0][2], 4.0, 2.0);
    assert_point(polylines[0][3], 2.0, 4.0);

    transform(&mut polylines, Scale::Same, Scale::Factor(0.5));
    assert_point(polylines[0][0], 0.0, 1.0);
    assert_point(polylines[0][1], 1.0, 0.0);
    assert_point(polylines[0][2], 2.0, 1.0);
    assert_point(polylines[0][3], 1.0, 2.0);

    transform(&mut polylines, Scale::Size(3.0), Scale::Factor(1.0));
    assert_point(polylines[0][0], 0.0, 1.0);
    assert_point(polylines[0][1], 1.5, 0.0);
    assert_point(polylines[0][2], 3.0, 1.0);
    assert_point(polylines[0][3], 1.5, 2.0);

    transform(&mut polylines, Scale::Factor(1.0), Scale::Size(5.0));
    assert_point(polylines[0][0], 0.0, 2.5);
    assert_point(polylines[0][1], 1.5, 0.0);
    assert_point(polylines[0][2], 3.0, 2.5);
    assert_point(polylines[0][3], 1.5, 5.0);
}

#[cfg(test)]
fn rotate_around(point: Point, center: Point, from_direction: Point, to_direction: Point) -> Point {
    let shifted = Point {
        x: point.x - center.x,
        y: point.y - center.y,
    };

    let cosine = from_direction.x * to_direction.x + from_direction.y * to_direction.y;
    let sine = from_direction.x * to_direction.y - from_direction.y * to_direction.x;

    let rotated = Point {
        x: shifted.x * cosine - shifted.y * sine,
        y: shifted.x * sine + shifted.y * cosine,
    };

    let unshifted = Point {
        x: rotated.x + center.x,
        y: rotated.y + center.y,
    };

    unshifted
}

#[cfg(test)]
#[test]
fn test_rotate_around() {
    let right = Point { x: 1.0, y: 0.0 };
    let up = Point { x: 0.0, y: 1.0 };
    let center = Point { x: 1.0, y: 1.0 };
    assert_point(
        rotate_around(Point { x: 1.0, y: 0.0 }, center, right, up),
        2.0,
        1.0,
    );
    assert_point(
        rotate_around(Point { x: 1.0, y: 2.0 }, center, up, right),
        2.0,
        1.0,
    );
}

#[cfg(test)]
fn arc_to(polyline: &mut Polyline, center: Point, to_point: Point, tolerance: f64) {
    let from_point = polyline.last().unwrap();
    let from_vector = Point {
        x: from_point.x - center.x,
        y: from_point.y - center.y,
    };

    let to_vector = Point {
        x: to_point.x - center.x,
        y: to_point.y - center.y,
    };

    let dot = from_vector.x * to_vector.x + from_vector.y * to_vector.y;
    let from_radius = (from_vector.x * from_vector.x + from_vector.y * from_vector.y).sqrt();
    let to_radius = (to_vector.x * to_vector.x + to_vector.y * to_vector.y).sqrt();
    let radius = (from_radius + to_radius) / 2.0;
    if radius < tolerance {
        return;
    }
    let det = from_vector.x * to_vector.y - from_vector.y * to_vector.x;

    let total_angle = det.atan2(dot);
    let max_step_angle = 2.0 * (1.0 - tolerance / radius).acos();
    let steps_amount = (total_angle.abs() / max_step_angle).ceil();
    let step_angle = total_angle / steps_amount;
    let steps_count = steps_amount as usize;

    for step_index in 1..steps_count {
        let angle = step_angle * step_index as f64;

        let cosine = angle.cos();
        let sine = angle.sin();

        let step_vector = Point {
            x: from_vector.x * cosine - from_vector.y * sine,
            y: from_vector.x * sine + from_vector.y * cosine,
        };

        let step_point = Point {
            x: center.x + step_vector.x,
            y: center.y + step_vector.y,
        };

        polyline.push(step_point);
    }

    polyline.push(to_point);
}

#[cfg(test)]
#[test]
fn test_arc_to() {
    let center = Point { x: 1.0, y: 1.0 };
    let bottom = Point { x: 1.0, y: 0.0 };
    let left = Point { x: 0.0, y: 1.0 };
    let right = Point { x: 2.0, y: 1.0 };
    let tolerance = 0.05;
    let small = 0.5;
    let large = 0.75_f64.sqrt();

    let mut left_polyline = vec![bottom];
    arc_to(&mut left_polyline, center, left, tolerance);
    assert!(left_polyline.len() == 4);
    assert_point(left_polyline[0], bottom.x, bottom.y);
    assert_point(left_polyline[1], 1.0 - small, 1.0 - large);
    assert_point(left_polyline[2], 1.0 - large, 1.0 - small);
    assert_point(left_polyline[3], left.x, left.y);

    let mut right_polyline = vec![bottom];
    arc_to(&mut right_polyline, center, right, tolerance);
    assert!(right_polyline.len() == 4);
    assert_point(right_polyline[0], bottom.x, bottom.y);
    assert_point(right_polyline[1], 1.0 + small, 1.0 - large);
    assert_point(right_polyline[2], 1.0 + large, 1.0 - small);
    assert_point(right_polyline[3], right.x, right.y);
}

/// The state we maintain for one of the paths.
#[derive(Debug)]
pub struct PathState<'a> {
    /// The current contact with the other path.
    pub contact_point: ContactPoint,

    /// A vector of the points along the path.
    pub points: &'a [Point],

    /// A vector of the directions of the path at each point.
    pub directions: &'a [Point],

    /// A vector of the lengths of the path segments at each point.
    pub lengths: &'a [f64],
}

/// Compute the path of the traced point as we rotate one path around another.
fn spiropath_polyline(
    _stationary_state: PathState,
    _rotating_state: PathState,
    _rotating_point: Point,
    _tolerance: f64,
) -> Polyline {
    _stationary_state.points.to_vec()
}

/// Generate spiropath(s) by rotating one shape around another.
pub fn spiropath(
    mut stationary: Polyline,
    mut rotating: Polyline,
    mut options: Options,
) -> Vec<Polyline> {
    if !is_clockwise(&stationary) {
        stationary.reverse();
    }

    match (options.location, is_clockwise(&rotating)) {
        (Location::Inside, false) => rotating.reverse(),
        (Location::Outside, true) => rotating.reverse(), // NOT TESTED
        _ => {}
    }

    let mut stationary_lengths = lengths(&stationary);
    let stationary_scale = options.stationary_teeth as f64 / stationary_lengths.iter().sum::<f64>();
    scale_lengths(&mut stationary_lengths, stationary_scale);
    scale_polyline(&mut stationary, stationary_scale, stationary_scale);
    let stationary_directions = directions(&stationary);

    let mut polylines: Vec<Polyline> = vec![];
    if options.include_stationary {
        polylines.push(stationary.clone());
    }

    let mut rotating_lengths = lengths(&rotating);
    let rotating_scale = options.rotating_teeth as f64 / rotating_lengths.iter().sum::<f64>();
    scale_lengths(&mut rotating_lengths, rotating_scale);
    scale_polyline(&mut rotating, rotating_scale, rotating_scale);
    scale_point(&mut options.traced_point, rotating_scale, rotating_scale);
    let rotating_directions = directions(&rotating);

    let rotating_top_most_index =
        find_top_most_point_index(&rotating, options.initial_rotating_angle);
    let rotating_contact_point = ContactPoint {
        index: rotating_top_most_index,
        offset: 0.0,
    };

    for initial_offset in options.initial_offsets.iter() {
        let stationary_state = PathState {
            contact_point: find_offset_contact_point(
                &stationary_lengths,
                *initial_offset,
                options.tolerance,
            ),
            points: &stationary,
            directions: &stationary_directions,
            lengths: &stationary_lengths,
        };

        let rotating_state = PathState {
            contact_point: rotating_contact_point,
            points: &rotating,
            directions: &rotating_directions,
            lengths: &rotating_lengths,
        };

        polylines.push(spiropath_polyline(
            stationary_state,
            rotating_state,
            options.traced_point,
            options.tolerance,
        ));
    }

    transform(&mut polylines, options.x_scale, options.y_scale);

    polylines
}

#[cfg(test)]
#[test]
fn test_spiropath() {
    let line = vec![Point { x: 0.0, y: 0.0 }, Point { x: 1.0, y: 0.0 }];

    let square = vec![
        Point { x: 0.0, y: 0.0 },
        Point { x: 1.0, y: 0.0 },
        Point { x: 1.0, y: 1.0 },
        Point { x: 0.0, y: 1.0 },
    ];

    let mut options = Options {
        location: Location::Inside,
        stationary_teeth: 4,
        rotating_teeth: 1,
        tolerance: 0.01,
        traced_point: Point { x: 0.0, y: 0.0 },
        initial_rotating_angle: 0.0,
        include_stationary: false,
        initial_offsets: vec![],
        x_scale: Scale::Same,
        y_scale: Scale::Same,
    };
    options.validate();

    let polylines = spiropath(square, line, options);

    let mut output = Vec::new();
    print_svg_polylines(&polylines, &mut output);

    let actual = from_utf8(&output).unwrap();
    let expected = "\
        <svg width='1pt', height='1pt', xmlns='http://www.w3.org/2000/svg'>\n\
        <path d='\n\
        M 0pt 1pt\n\
        L 1pt 1pt\n\
        L 1pt 0pt\n\
        L 0pt 0pt\n\
        Z\n\
        '/>\n\
        </svg>\n\
    ";

    assert!(
        actual == expected,
        "The actual results:\n\
             ---\n\
             {}\
             ...\n\
             Are different from the expected results:\n\
             ---\n\
             {}\
             ...",
        actual,
        expected
    );
}
