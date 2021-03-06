use ordered_float::OrderedFloat;
use std::cmp::min;
use std::cmp::Ordering;
use std::f64::consts::PI;
use std::fs::read_to_string;
use std::io::Write;
use std::mem::swap;
use svg2polylines::parse as parse_svg;
use svg2polylines::CoordinatePair as Point;
use svg2polylines::Polyline;

/// Return the distance between two points.
pub fn points_distance(left_point: Point, right_point: Point) -> f64 {
    let delta_x = left_point.x - right_point.x;
    let delta_y = left_point.y - right_point.y;
    delta_x.hypot(delta_y)
}

/// Measure the lengths of each of the segments of the polyline.
pub fn polyline_lengths(polyline: &[Point]) -> Vec<f64> {
    (0..polyline.len())
        .map(|prev_index| {
            let next_index = (prev_index + 1) % polyline.len();
            points_distance(polyline[prev_index], polyline[next_index])
        })
        .collect()
}

/// Scale a polyline by a factor.
pub fn scale_lengths(lengths: &mut Vec<f64>, factor: f64) {
    for length in lengths.iter_mut() {
        *length *= factor;
    }
}

/// Round lengths to integers.
pub fn round_lengths(lengths: &[f64], factor: f64) -> Vec<usize> {
    lengths
        .iter()
        .map(|length| (length * factor).round() as usize)
        .collect()
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
    assert!(polyline_lengths(&square).iter().sum::<f64>() == 4.0);
    scale_polyline(&mut square, 2.0, 3.0);
    assert!(polyline_lengths(&square).iter().sum::<f64>() == 10.0);
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

/// Return the direction of the polyline at each of its points.
pub fn polyline_directions(polyline: &[Point]) -> Polyline {
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
    let triangle_directions = polyline_directions(&triangle);

    assert_point(triangle_directions[0], 1.0, 0.0);
    assert_point(triangle_directions[1], 0.0, 1.0);
    assert_point(triangle_directions[2], -0.5_f64.sqrt(), -0.5_f64.sqrt());
}

/// How to interpret 180 degree bends.
#[derive(Clone, Copy, Debug)]
pub enum Orientation {
    /// Assume clockwise.
    Clockwise,

    /// Assume Widdershins.
    Widdershins,
}

/// Orient a 180 degree turn.
pub fn orient_angle(angle: f64, orientation: Orientation) -> f64 {
    match orientation {
        Orientation::Clockwise if angle < -PI + 1e-4 => angle + 2.0 * PI,
        Orientation::Widdershins if angle > PI - 1e-4 => angle - 2.0 * PI,
        _ => angle,
    }
}

/// Return the bending angle at each point along a polyline.
pub fn directions_bends(directions: &[Point], orientation: Orientation) -> Vec<f64> {
    let mut result = vec![0.0; directions.len()];
    for next_index in 0..directions.len() {
        let prev_index = (next_index + directions.len() - 1) % directions.len();
        let next_direction = directions[next_index];
        let prev_direction = directions[prev_index];
        let cosine = next_direction.x * prev_direction.x + next_direction.y * prev_direction.y;
        let sine = next_direction.x * prev_direction.y - next_direction.y * prev_direction.x;
        let angle = sine.atan2(cosine);
        result[next_index] = orient_angle(angle, orientation);
    }
    result
}

#[cfg(test)]
#[test]
fn test_bends() {
    let mut triangle = vec![
        Point { x: 0.0, y: 0.0 },
        Point { x: 1.0, y: 1.0 },
        Point { x: 1.0, y: 0.0 },
    ];
    let mut triangle_directions = polyline_directions(&triangle);
    let mut triangle_bends = directions_bends(&triangle_directions, Orientation::Clockwise);

    assert_float_absolute_eq!(triangle_bends[0] * 180.0 / PI, 135.0, 1e-6);
    assert_float_absolute_eq!(triangle_bends[1] * 180.0 / PI, 135.0, 1e-6);
    assert_float_absolute_eq!(triangle_bends[2] * 180.0 / PI, 90.0, 1e-6);

    triangle.reverse();
    triangle_directions = polyline_directions(&triangle);
    triangle_bends = directions_bends(&triangle_directions, Orientation::Clockwise);

    assert_float_absolute_eq!(triangle_bends[0] * 180.0 / PI, -90.0, 1e-6);
    assert_float_absolute_eq!(triangle_bends[1] * 180.0 / PI, -135.0, 1e-6);
    assert_float_absolute_eq!(triangle_bends[2] * 180.0 / PI, -135.0, 1e-6);

    let mut house = vec![
        Point { x: 0.0, y: 0.0 },
        Point { x: 0.0, y: 1.0 },
        Point { x: 1.0, y: 2.0 },
        Point { x: 2.0, y: 1.0 },
        Point { x: 2.0, y: 0.0 },
    ];
    let mut house_directions = polyline_directions(&house);
    let mut house_bends = directions_bends(&house_directions, Orientation::Clockwise);

    assert_float_absolute_eq!(house_bends[0] * 180.0 / PI, 90.0, 1e-6);
    assert_float_absolute_eq!(house_bends[1] * 180.0 / PI, 45.0, 1e-6);
    assert_float_absolute_eq!(house_bends[2] * 180.0 / PI, 90.0, 1e-6);
    assert_float_absolute_eq!(house_bends[3] * 180.0 / PI, 45.0, 1e-6);
    assert_float_absolute_eq!(house_bends[4] * 180.0 / PI, 90.0, 1e-6);

    house.reverse();
    house_directions = polyline_directions(&house);
    house_bends = directions_bends(&house_directions, Orientation::Clockwise);

    assert_float_absolute_eq!(house_bends[0] * 180.0 / PI, -90.0, 1e-6);
    assert_float_absolute_eq!(house_bends[1] * 180.0 / PI, -45.0, 1e-6);
    assert_float_absolute_eq!(house_bends[2] * 180.0 / PI, -90.0, 1e-6);
    assert_float_absolute_eq!(house_bends[3] * 180.0 / PI, -45.0, 1e-6);
    assert_float_absolute_eq!(house_bends[4] * 180.0 / PI, -90.0, 1e-6);

    let line = vec![Point { x: 0.0, y: 0.0 }, Point { x: 1.0, y: 0.0 }];

    let line_directions = polyline_directions(&line);
    let line_bends = directions_bends(&line_directions, Orientation::Clockwise);
    assert_float_absolute_eq!(line_bends[0] * 180.0 / PI, 180.0, 1e-6);
    assert_float_absolute_eq!(line_bends[1] * 180.0 / PI, 180.0, 1e-6);
}

/// Test whether a pre-processed polyline has a clockwise orientation.
pub fn is_clockwise(polyline: &[Point]) -> bool {
    let minimal_y = *polyline
        .iter()
        .map(|point| OrderedFloat(point.y))
        .min()
        .unwrap();
    let mut sum = 0.0;
    for prev_index in 0..polyline.len() {
        let next_index = (prev_index + 1) % polyline.len();
        let prev_point = polyline[prev_index];
        let next_point = polyline[next_index];
        sum +=
            (next_point.x - prev_point.x) * (next_point.y + prev_point.y - 2.0 * minimal_y + 100.0);
    }
    sum >= 0.0
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

/// Return the minimal coordinates in some polylines.
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

/// Return the maximal coordinates in some polylines.
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

/// Scale and move the points to fit in a (0,0) -> (w,h) bounding box.
pub fn transform(polylines: &mut [Polyline], x_scale: Scale, y_scale: Scale) {
    let minimal_point = minimal_coordinates(polylines);
    let maximal_point = maximal_coordinates(polylines);
    let size = Point {
        x: maximal_point.x - minimal_point.x,
        y: maximal_point.y - minimal_point.y,
    };

    let maybe_x_factor = match x_scale {
        Scale::Size(x_size) => Some(if size.x > 0.0 { x_size / size.x } else { 1.0 }),
        Scale::Factor(factor) => Some(factor),
        Scale::Same => None,
    };

    let maybe_y_factor = match y_scale {
        Scale::Size(y_size) => Some(if size.y > 0.0 { y_size / size.y } else { 1.0 }),
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

pub fn prune_repeated_points(polyline: &mut Polyline, tolerance: f64) {
    let mut kept_count = 1;
    for next_index in 1..polyline.len() {
        if points_distance(polyline[kept_count - 1], polyline[next_index]) <= tolerance {
            continue;
        }
        if kept_count < next_index {
            polyline[kept_count] = polyline[next_index];
            kept_count += 1;
        }
    }

    if kept_count > 1 && points_distance(polyline[0], polyline[kept_count - 1]) <= tolerance {
        kept_count -= 1;
    }

    polyline.truncate(kept_count);
}

/// Load a polyline from a file and convert it to a polyline.
pub fn load_polyline_from_svg_file(path: &str, tolerance: f64) -> Polyline {
    let string = read_to_string(path).unwrap_or_else(|_| panic!("reading file: {}", path));

    let mut polylines = parse_svg(&string).unwrap_or_else(|_| panic!("parsing file: {}", path));

    if polylines.is_empty() {
        panic!("no SVG paths in file: {}", path); // NOT TESTED
    }

    if polylines.len() > 1 {
        panic!("too many SVG paths in file: {}", path); // NOT TESTED
    }

    let mut polyline = polylines.pop().unwrap();

    if polyline.len() < 2 {
        panic!("too few points in SVG path in file: {}", path); // NOT TESTED
    }

    prune_repeated_points(&mut polyline, tolerance);

    if polyline.len() == 1 {
        panic!("zero length path in file: {}", path); // NOT TESTED
    }

    polyline
}

fn print_svg_polyline(polyline: &[Point], id: &str, width: f64, output: &mut dyn Write) {
    writeln!(
        output,
        "<path id='{}' fill='none' stroke='black' stroke-width='{}' d='",
        id, width
    )
    .unwrap();
    let mut command = "M";
    for point in polyline {
        writeln!(output, "{} {} {}", command, point.x, point.y).unwrap();
        command = "L";
    }
    writeln!(output, "Z").unwrap();
    writeln!(output, "'/>").unwrap();
}

/// Print a vector of polylines as an SVG file.
pub fn print_svg_polylines(polylines: &[Polyline], ids: &[String], output: &mut dyn Write) {
    let maximal_point = maximal_coordinates(polylines);
    writeln!(
        output,
        "<svg width='{}pt' height='{}pt' xmlns='http://www.w3.org/2000/svg'>",
        maximal_point.x, maximal_point.y
    )
    .unwrap();

    writeln!(output, "<g transform='scale(1.333333 1.333333)'>").unwrap();

    let width = ((maximal_point.x / 1000.0) * (maximal_point.y / 1000.0)).sqrt();
    for (polyline, id) in polylines.iter().zip(ids.iter()) {
        print_svg_polyline(polyline, id, width, output);
    }

    writeln!(output, "</g>").unwrap();
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

    /// Whether to mirror the rotating shape.
    pub mirror_rotating: bool,

    /// Whether to include the stationary shape in the output.
    pub include_stationary: bool,

    /// The initial offsets of the rotating shape.
    pub initial_offsets: Vec<f64>,

    /// How to scale the X axis of the output.
    pub x_scale: Scale,

    /// How to scale the Y axis of the output.
    pub y_scale: Scale,

    /// Duration of animation, or 0 for nonw.
    pub duration: f64,
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
            panic!("no offsets given"); // NOT TESTED
        }
    }
}

/// A point along a polyline where it is in contact with another.
#[derive(Clone, Copy, Debug)]
pub struct ContactPoint {
    /// The index of the point in the polyline.
    pub index: usize,

    /// How far to move along the line segment starting at that point.
    pub offset: usize,
}

/// Find a point along a polyline by its offset (up to some tolerance).
fn find_offset_contact_point(lengths: &[usize], mut offset: usize) -> ContactPoint {
    let mut index = 0;

    while offset > lengths[index] {
        offset -= lengths[index];
        index = (index + 1) % lengths.len();
    }

    ContactPoint { index, offset }
}

fn traced_position(
    point: Point,
    from_center: Point,
    to_center: Point,
    from_direction: Point,
    to_direction: Point,
) -> Point {
    let shifted = Point {
        x: point.x - from_center.x,
        y: point.y - from_center.y,
    };

    let cosine = from_direction.x * to_direction.x + from_direction.y * to_direction.y;
    let mut sine = from_direction.x * to_direction.y - from_direction.y * to_direction.x;

    let angle = sine.atan2(cosine);
    if angle < -PI + 1e-4 {
        sine *= -1.0;
    }

    let rotated = Point {
        x: shifted.x * cosine - shifted.y * sine,
        y: shifted.x * sine + shifted.y * cosine,
    };

    Point {
        x: rotated.x + to_center.x,
        y: rotated.y + to_center.y,
    }
}

#[cfg(test)]
#[test]
fn test_traced_position() {
    let right = Point { x: 1.0, y: 0.0 };
    let up = Point { x: 0.0, y: 1.0 };
    let from_center = Point { x: 1.0, y: 1.0 };
    let to_center = Point { x: 11.0, y: 11.0 };
    assert_point(
        traced_position(Point { x: 1.0, y: 0.0 }, from_center, to_center, right, up),
        12.0,
        11.0,
    );
    assert_point(
        traced_position(Point { x: 1.0, y: 2.0 }, from_center, to_center, up, right),
        12.0,
        11.0,
    );
}

fn arc_to(
    polyline: &mut Polyline,
    center: Point,
    from_point: Point,
    to_point: Point,
    arc_angle: f64,
    tolerance: f64,
) {
    let from_vector = Point {
        x: from_point.x - center.x,
        y: from_point.y - center.y,
    };

    let to_vector = Point {
        x: to_point.x - center.x,
        y: to_point.y - center.y,
    };

    let from_radius = (from_vector.x * from_vector.x + from_vector.y * from_vector.y).sqrt();
    let to_radius = (to_vector.x * to_vector.x + to_vector.y * to_vector.y).sqrt();
    let radius = (from_radius + to_radius) / 2.0;
    if radius < tolerance {
        return;
    }

    let max_step_angle = 2.0 * (1.0 - tolerance / radius).acos();
    let steps_amount = (arc_angle.abs() / max_step_angle).ceil();
    let step_angle = arc_angle / steps_amount;
    let steps_count = steps_amount as usize;

    polyline.push(from_point);
    for step_index in 1..steps_count {
        let angle = step_angle * step_index as f64;

        let cosine = angle.cos();
        let sine = angle.sin();

        let step_vector = Point {
            x: from_vector.y * sine + from_vector.x * cosine,
            y: from_vector.y * cosine - from_vector.x * sine,
        };

        let step_point = Point {
            x: center.x + step_vector.x,
            y: center.y + step_vector.y,
        };

        polyline.push(step_point);
    }
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

    let mut left_polyline: Polyline = vec![];
    arc_to(
        &mut left_polyline,
        center,
        bottom,
        left,
        90.0 * PI / 180.0,
        tolerance,
    );
    assert!(left_polyline.len() == 3);
    assert_point(left_polyline[0], bottom.x, bottom.y);
    assert_point(left_polyline[1], 1.0 - small, 1.0 - large);
    assert_point(left_polyline[2], 1.0 - large, 1.0 - small);

    let mut right_polyline: Polyline = vec![];
    arc_to(
        &mut right_polyline,
        center,
        bottom,
        right,
        -90.0 * PI / 180.0,
        tolerance,
    );
    assert!(right_polyline.len() == 3);
    assert_point(right_polyline[0], bottom.x, bottom.y);
    assert_point(right_polyline[1], 1.0 + small, 1.0 - large);
    assert_point(right_polyline[2], 1.0 + large, 1.0 - small);
}

/// The state we maintain for one of the paths.
#[derive(Clone, Copy, Debug)]
pub struct PathState<'a> {
    /// The current contact with the other polyline.
    pub contact_point: ContactPoint,

    /// A vector of the points along the polyline.
    pub points: &'a [Point],

    /// A vector of the directions of the polyline at each point.
    pub directions: &'a [Point],

    /// A vector of the lengths of the polyline segments at each point.
    pub lengths: &'a [usize],

    /// A vector of the bending angles of the polyline segments at each point.
    pub bends: &'a [f64],
}

fn state_center_and_position(
    stationary_state: &PathState,
    rotated_state: &PathState,
    traced_point: Point,
) -> (Point, Point) {
    let from_direction = rotated_state.directions[rotated_state.contact_point.index];
    let to_direction = stationary_state.directions[stationary_state.contact_point.index];

    let from_center = Point {
        x: rotated_state.points[rotated_state.contact_point.index].x
            + from_direction.x * rotated_state.contact_point.offset as f64,
        y: rotated_state.points[rotated_state.contact_point.index].y
            + from_direction.y * rotated_state.contact_point.offset as f64,
    };
    let to_center = Point {
        x: stationary_state.points[stationary_state.contact_point.index].x
            + to_direction.x * stationary_state.contact_point.offset as f64,
        y: stationary_state.points[stationary_state.contact_point.index].y
            + to_direction.y * stationary_state.contact_point.offset as f64,
    };

    let position = traced_position(
        traced_point,
        from_center,
        to_center,
        from_direction,
        to_direction,
    );

    (to_center, position)
}

fn next_contact_point(
    stationary_state: &mut PathState,
    rotating_state: &mut PathState,
) -> (usize, f64) {
    let stationary_distance = stationary_state.lengths[stationary_state.contact_point.index]
        - stationary_state.contact_point.offset;
    let rotating_distance = rotating_state.lengths[rotating_state.contact_point.index]
        - rotating_state.contact_point.offset;

    match stationary_distance.cmp(&rotating_distance) {
        Ordering::Less => {
            stationary_state.contact_point.index =
                (stationary_state.contact_point.index + 1) % stationary_state.points.len();
            stationary_state.contact_point.offset = 0;

            rotating_state.contact_point.offset += stationary_distance;

            (
                stationary_distance,
                stationary_state.bends[stationary_state.contact_point.index],
            )
        }

        Ordering::Greater => {
            rotating_state.contact_point.index =
                (rotating_state.contact_point.index + 1) % rotating_state.points.len();
            rotating_state.contact_point.offset = 0;

            stationary_state.contact_point.offset += rotating_distance;

            (
                rotating_distance,
                -rotating_state.bends[rotating_state.contact_point.index],
            )
        }

        Ordering::Equal => {
            stationary_state.contact_point.index =
                (stationary_state.contact_point.index + 1) % stationary_state.points.len();
            stationary_state.contact_point.offset = 0;

            rotating_state.contact_point.index =
                (rotating_state.contact_point.index + 1) % rotating_state.points.len();
            rotating_state.contact_point.offset = 0;

            (
                stationary_distance,
                stationary_state.bends[stationary_state.contact_point.index]
                    - rotating_state.bends[rotating_state.contact_point.index],
            )
        }
    }
}

/// Compute the polyline of the traced point as we rotate one polyline around another.
pub fn spiropath_polyline(
    stationary_state: &PathState,
    rotating_state: &PathState,
    traced_point: Point,
    total_length: usize,
) -> Polyline {
    let mut current_stationary_state = *stationary_state;
    let mut current_rotating_state = *rotating_state;
    let mut polyline: Polyline = Vec::new();
    let mut contact_length: usize = 0;

    let (_, mut current_traced_point) = state_center_and_position(
        &current_stationary_state,
        &current_rotating_state,
        traced_point,
    );

    while contact_length < total_length {
        let (offset, angle) =
            next_contact_point(&mut current_stationary_state, &mut current_rotating_state);
        contact_length += offset;

        let (next_center, next_traced_point) = state_center_and_position(
            &current_stationary_state,
            &current_rotating_state,
            traced_point,
        );

        arc_to(
            &mut polyline,
            next_center,
            current_traced_point,
            next_traced_point,
            angle,
            1.0,
        );

        current_traced_point = next_traced_point;
    }

    assert!(contact_length == total_length);
    polyline
}

fn gcd(mut n: usize, mut m: usize) -> usize {
    if n > m {
        swap(&mut n, &mut m);
    }
    while n > 0 && m > 0 {
        m %= n;
        swap(&mut n, &mut m);
    }
    n + m
}

fn lcm(n: usize, m: usize) -> usize {
    n * m / gcd(n, m)
}

#[test]
fn test_lcm() {
    assert_eq!(lcm(8, 12), 24);
}

fn round_lengths_scale(
    stationary_raw_lengths: &[f64],
    stationary_teeth: usize,
    rotating_raw_lengths: &[f64],
    rotating_teeth: usize,
    mut tolerance: f64,
) -> (f64, usize) {
    let stationary_min_length: f64 = *stationary_raw_lengths
        .iter()
        .map(|length| OrderedFloat(*length))
        .min()
        .unwrap();

    let rotating_min_length: f64 = *rotating_raw_lengths
        .iter()
        .map(|length| OrderedFloat(*length))
        .min()
        .unwrap();

    tolerance = *min(OrderedFloat(tolerance), OrderedFloat(stationary_min_length));
    tolerance = *min(OrderedFloat(tolerance), OrderedFloat(rotating_min_length));

    let mut scale = (10.0 / tolerance).ceil();
    let factor = 1.0001;
    let mut attempts = 1;
    loop {
        let stationary_sum = stationary_raw_lengths
            .iter()
            .map(|length| (length * scale).round() as usize)
            .sum::<usize>();

        let rotating_sum = rotating_raw_lengths
            .iter()
            .map(|length| (length * scale).round() as usize)
            .sum::<usize>();

        let stationary_remainder = stationary_sum % stationary_teeth;
        let rotating_remainder = rotating_sum % rotating_teeth;
        let stationary_scale = stationary_sum / stationary_teeth;
        let rotating_scale = rotating_sum / rotating_teeth;

        if stationary_remainder == 0
            && rotating_remainder == 0
            && stationary_scale == rotating_scale
        {
            return (
                scale,
                lcm(stationary_teeth, rotating_teeth) * rotating_scale,
            );
        }

        scale *= factor;

        attempts += 1;
        if attempts > 10000 {
            panic!("failed to scale lengths to integer values"); // NOT TESTED
        }
    }
}

/// Generate spiropath(s) by rotating one shape around another.
pub fn spiropath(
    mut stationary: Polyline,
    mut rotating: Polyline,
    mut options: Options,
) -> (Vec<Polyline>, Vec<String>) {
    if !is_clockwise(&stationary) {
        stationary.reverse();
    }

    if options.mirror_rotating {
        // BEGIN NOT TESTED
        options.traced_point.x *= -1.0;
        for point in rotating.iter_mut() {
            point.x *= -1.0;
        }
        // END NOT TESTED
    }

    match (options.location, is_clockwise(&rotating)) {
        (Location::Inside, false) => {
            rotating.reverse(); // NOT TESTED
        }
        (Location::Outside, true) => {
            rotating.reverse();
        }
        _ => {}
    };

    let mut stationary_raw_lengths = polyline_lengths(&stationary);
    let mut rotating_raw_lengths = polyline_lengths(&rotating);

    let rotating_directions = polyline_directions(&rotating);
    let stationary_directions = polyline_directions(&stationary);

    let rotating_orientation = match options.location {
        Location::Inside => Orientation::Clockwise,
        Location::Outside => Orientation::Widdershins,
    };

    let rotating_bends = directions_bends(&rotating_directions, rotating_orientation);
    let stationary_bends = directions_bends(&stationary_directions, Orientation::Clockwise);

    let stationary_scale =
        options.stationary_teeth as f64 / stationary_raw_lengths.iter().sum::<f64>();
    let rotating_scale = options.rotating_teeth as f64 / rotating_raw_lengths.iter().sum::<f64>();

    scale_lengths(&mut stationary_raw_lengths, stationary_scale);
    scale_lengths(&mut rotating_raw_lengths, rotating_scale);

    let (round_scale, total_length) = round_lengths_scale(
        &stationary_raw_lengths,
        options.stationary_teeth,
        &rotating_raw_lengths,
        options.rotating_teeth,
        options.tolerance,
    );

    let stationary_lengths = round_lengths(&stationary_raw_lengths, round_scale);
    let rotating_lengths = round_lengths(&rotating_raw_lengths, round_scale);

    let mut ids: Vec<String> = Vec::new();
    let mut polylines: Vec<Polyline> = vec![];
    if options.include_stationary {
        let mut polyline = stationary.clone();
        scale_polyline(&mut polyline, stationary_scale, stationary_scale);
        polylines.push(polyline);
        ids.push("stationary".to_owned());
    }

    scale_polyline(
        &mut stationary,
        round_scale * stationary_scale,
        round_scale * stationary_scale,
    );
    scale_polyline(
        &mut rotating,
        round_scale * rotating_scale,
        round_scale * rotating_scale,
    );
    scale_point(
        &mut options.traced_point,
        round_scale * rotating_scale,
        round_scale * rotating_scale,
    );

    let rotating_top_most_index =
        find_top_most_point_index(&rotating, options.initial_rotating_angle);
    let rotating_contact_point = ContactPoint {
        index: rotating_top_most_index,
        offset: 0,
    };

    for initial_offset in options.initial_offsets.iter() {
        let stationary_state = PathState {
            contact_point: find_offset_contact_point(
                &stationary_lengths,
                (*initial_offset * round_scale).round() as usize,
            ),
            points: &stationary,
            directions: &stationary_directions,
            lengths: &stationary_lengths,
            bends: &stationary_bends,
        };

        let rotating_state = PathState {
            contact_point: rotating_contact_point,
            points: &rotating,
            directions: &rotating_directions,
            lengths: &rotating_lengths,
            bends: &rotating_bends,
        };

        ids.push(format!("{:?}-{}", options.location, initial_offset));
        let mut polyline = spiropath_polyline(
            &stationary_state,
            &rotating_state,
            options.traced_point,
            total_length,
        );
        scale_polyline(&mut polyline, 1.0 / round_scale, 1.0 / round_scale);
        polylines.push(polyline);
    }

    transform(&mut polylines, options.x_scale, options.y_scale);

    (polylines, ids)
}
