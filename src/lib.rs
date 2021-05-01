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

use ordered_float::OrderedFloat;
use std::cmp::Ordering;
use std::fs::read_to_string;
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
    let mut polyline = vec![
        Point { x: 0.0, y: 0.0 },
        Point { x: 1.0, y: 0.0 },
        Point { x: 1.0, y: 1.0 },
    ];
    assert!(!is_clockwise(&polyline));
    polyline.reverse();
    assert!(is_clockwise(&polyline));
}

/// Measure the lengths of each of the segments of the polyline.
pub fn lengths(polyline: &[Point]) -> Vec<f64> {
    (0..polyline.len())
        .map(|prev_index| {
            let next_index = (prev_index + 1) % polyline.len();
            let prev_point = polyline[prev_index];
            let next_point = polyline[next_index];
            let delta_x = next_point.x - prev_point.x;
            let delta_y = next_point.y - prev_point.y;
            delta_x.hypot(delta_y)
        })
        .collect()
}

/// Scale a polyline by a factor.
pub fn scale_lengths(lengths: &mut Vec<f64>, factor: f64) {
    for length in lengths.iter_mut() {
        *length *= factor;
    }
}

/// Scale a polyline by a factor.
pub fn scale_polyline(polyline: &mut Polyline, x_factor: f64, y_factor: f64) {
    for point in polyline.iter_mut() {
        point.x *= x_factor;
        point.y *= y_factor;
    }
}

#[cfg(test)]
#[test]
fn test_scale() {
    let mut polyline = vec![
        Point { x: 0.0, y: 0.0 },
        Point { x: 1.0, y: 0.0 },
        Point { x: 1.0, y: 1.0 },
        Point { x: 0.0, y: 1.0 },
    ];
    assert!(lengths(&polyline).iter().sum::<f64>() == 4.0);
    scale_polyline(&mut polyline, 2.0, 3.0);
    assert!(lengths(&polyline).iter().sum::<f64>() == 10.0);
}

/// Load a path from a file and convert it to a polyline.
pub fn load_polyline_from_svg_file(path: &str) -> Polyline {
    let string = read_to_string(path).unwrap_or_else(|_| panic!("reading file: {}", path));

    let mut polylines = parse_svg(&string).unwrap_or_else(|_| panic!("parsing file: {}", path));

    if polylines.is_empty() {
        panic!("no SVG paths in file: {}", path);
    }

    if polylines.len() > 1 {
        panic!("too many SVG paths in file: {}", path);
    }

    let polyline = polylines.pop().unwrap();

    if polyline.len() < 2 {
        panic!("too few points in SVG path in file: {}", path);
    }

    if lengths(&polyline).iter().sum::<f64>() == 0.0 {
        panic!("zero length path in file: {}", path);
    }

    polyline
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
            panic!("stationary-teeth is zero");
        }

        if self.rotating_teeth == 0 {
            panic!("rotating-teeth is zero");
        }

        if self.tolerance <= 0.0 {
            panic!("tolerance {} is not positive", self.tolerance);
        }

        if self.initial_offsets.is_empty() {
            self.initial_offsets.push(0.0);
        }
    }
}

/// Generate spiropath(s) by rotating one shape around another.
pub fn spiropath(
    mut stationary: Polyline,
    mut rotating: Polyline,
    options: &Options,
) -> Vec<Polyline> {
    if !is_clockwise(&stationary) {
        stationary.reverse();
    }

    match (options.location, is_clockwise(&rotating)) {
        (Location::Outside, false) => rotating.reverse(),
        (Location::Inside, true) => rotating.reverse(),
        _ => {}
    }

    let mut stationary_lengths = lengths(&stationary);
    let stationary_scale = options.stationary_teeth as f64 / stationary_lengths.iter().sum::<f64>();
    scale_lengths(&mut stationary_lengths, stationary_scale);
    scale_polyline(&mut stationary, stationary_scale, stationary_scale);

    let mut rotating_lengths = lengths(&rotating);
    let rotating_scale = options.rotating_teeth as f64 / rotating_lengths.iter().sum::<f64>();
    scale_lengths(&mut rotating_lengths, rotating_scale);
    scale_polyline(&mut rotating, rotating_scale, rotating_scale);

    let mut polylines: Vec<Polyline> = Vec::new();
    if options.include_stationary {
        polylines.push(stationary);
    }

    panic!("not implemented");

    // TODOX polylines
}
