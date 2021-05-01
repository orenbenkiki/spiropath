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

use clap::crate_version;
use clap::App;
use clap::Arg;
use clap::ArgMatches;
use spiropath::load_polyline_from_svg_file;
use spiropath::Location;
use spiropath::Options;
use spiropath::Scale;
use std::f64::consts::PI;
use svg2polylines::parse as parse_svg;
use svg2polylines::CoordinatePair as Point;
use svg2polylines::Polyline;

fn main() {
    let arg_matches = app().get_matches();
    let mut options = arg_options(&arg_matches);
    options.validate();
    let stationary = parse_polyline(&arg_matches, "stationary");
    let rotating = parse_polyline(&arg_matches, "rotating");
    let output_path = arg_matches.value_of("output");
    eprintln!("Options: {:?}", options);
    eprintln!("Stationary: {:?}", stationary);
    eprintln!("Rotating: {:?}", rotating);
    eprintln!("Output: {:?}", output_path);
}

fn app() -> App<'static, 'static> {
    App::new("spiropath")
        .about("\nGeneralized Spirograph using arbitrary paths.")
        .version(crate_version!())
        .arg(
            Arg::with_name("stationary")
                .long("stationary")
                .short("s")
                .value_name("FILE or COUNT")
                .help("SVG file containing the stationary path, or the number of sides of a radius-1 polygon \
                      (1 - circle, 2 - line, 3 - triangle, etc.)")
                .required(true)
        )
        .arg(
            Arg::with_name("rotating")
                .long("rotating")
                .short("r")
                .value_name("FILE or COUNT")
                .help("SVG file containing the rotating path, or the number of sides of a radius-1 polygon \
                      (1 - circle, 2 - line, 3 - triangle, etc.)")
                .required(true)
        )
        .arg(
            Arg::with_name("output")
                .long("output")
                .short("o")
                .value_name("FILE")
                .help("SVG file to write the spiropath into")
        )
        .arg(
            Arg::with_name("location")
                .long("location")
                .short("L")
                .value_name("SIDE")
                .possible_values(&["inside", "outside"])
                .help("How to place the rotating path relative to the stationary path")
                .default_value("outside")
        )
        .arg(
            Arg::with_name("stationary-teeth")
                .long("stationary-teeth")
                .short("S")
                .value_name("COUNT")
                .help("Number of teeth on the stationary path")
                .default_value("47")
        )
        .arg(
            Arg::with_name("rotating-teeth")
                .long("rotating-teeth")
                .short("R")
                .value_name("COUNT")
                .help("Number of teeth on the rotating path")
                .default_value("7")
        )
        .arg(
            Arg::with_name("tolerance")
                .long("tolerance")
                .short("T")
                .value_name("FRACTION")
                .help("Linear path approximation tolerance (fraction of the teeth size)")
                .default_value("0.1")
        )
        .arg(
            Arg::with_name("point")
                .long("point")
                .short("P")
                .help("Coordinates of the traced point (relative to the rotating path)")
                .value_name("X,Y")
                .number_of_values(2)
                .require_delimiter(true)
                .default_value("0,0")
        )
        .arg(
            Arg::with_name("angle")
                .long("angle")
                .short("A")
                .value_name("DEGREES")
                .help("Initial angle (in degrees) of the rotating path")
                .default_value("0")
        )
        .arg(
            Arg::with_name("offset")
                .long("offset")
                .short("O")
                .value_name("FRACTION(S)")
                .help("Offset(s) of start position(s) of rotating path \
                      relative to top point of stationary path (fraction of the teeth size)")
                .multiple(true)
                .use_delimiter(true)
                .min_values(1)
                .default_value("0")
        )
        .arg(
            Arg::with_name("x-scale")
                .long("x-scale")
                .short("X")
                .value_name("SCALE")
                .help("Scaling of the output SVG (\"<size>pt\" / \"x<factor>\" / \"same\" as y-scale)")
                .default_value("x1")
        )
        .arg(
            Arg::with_name("y-scale")
                .long("y-scale")
                .short("Y")
                .value_name("SCALE")
                .help("Scaling of the output SVG (\"<size>pt\" / \"x<factor>\" / \"same\" as x-scale)")
                .default_value("x1")
        )
        .arg(
            Arg::with_name("include-stationary")
                .long("include-stationary")
                .short("I")
                .help("Include the (scaled) stationary path in the output")
        )
}

fn arg_options(arg_matches: &ArgMatches) -> Options {
    let location = parse_location(arg_matches, "location");
    let stationary_teeth = parse_count(arg_matches, "stationary-teeth");
    let rotating_teeth = parse_count(arg_matches, "rotating-teeth");
    let tolerance = parse_fraction(arg_matches, "tolerance");
    let traced_point = parse_point(arg_matches, "point");
    let initial_rotating_angle = parse_value(arg_matches, "angle");
    let include_stationary = arg_matches.is_present("include-stationary");
    let initial_offsets = parse_offsets(arg_matches, "offset", stationary_teeth as f64);
    let x_scale = parse_scale(arg_matches, "x-scale");
    let y_scale = parse_scale(arg_matches, "x-scale");
    Options {
        location,
        stationary_teeth,
        rotating_teeth,
        tolerance,
        traced_point,
        initial_rotating_angle,
        include_stationary,
        initial_offsets,
        x_scale,
        y_scale,
    }
}

fn parse_location(arg_matches: &ArgMatches, name: &str) -> Location {
    match arg_matches.value_of(name).unwrap() {
        "inside" => Location::Inside,
        "outside" => Location::Outside,
        _ => unreachable!(),
    }
}

fn parse_count(arg_matches: &ArgMatches, name: &str) -> usize {
    let value = arg_matches
        .value_of(name)
        .unwrap()
        .parse::<usize>()
        .unwrap_or_else(|error| {
            panic!(
                "{} in {}: {}",
                error,
                name,
                arg_matches.value_of(name).unwrap()
            )
        });
    if value == 0 {
        panic!("{} is zero", name);
    }
    value
}

fn parse_fraction(arg_matches: &ArgMatches, name: &str) -> f64 {
    let value = parse_value(arg_matches, name);
    if value <= 0.0 {
        panic!("{}: {} is not positive", name, value);
    }
    value
}

fn parse_value(arg_matches: &ArgMatches, name: &str) -> f64 {
    arg_matches
        .value_of(name)
        .unwrap()
        .parse::<f64>()
        .unwrap_or_else(|error| {
            panic!(
                "{} in {}: {}",
                error,
                name,
                arg_matches.value_of(name).unwrap()
            )
        })
}

fn parse_point(arg_matches: &ArgMatches, name: &str) -> Point {
    let coordinates: Vec<f64> = parse_values(arg_matches, name);
    assert!(coordinates.len() == 2);
    Point {
        x: coordinates[0],
        y: coordinates[1],
    }
}

fn parse_values(arg_matches: &ArgMatches, name: &str) -> Vec<f64> {
    arg_matches
        .values_of(name)
        .unwrap()
        .map(|string| {
            string
                .parse::<f64>()
                .unwrap_or_else(|error| panic!("{} in {}: {}", error, name, string))
        })
        .collect()
}

fn parse_offsets(arg_matches: &ArgMatches, name: &str, stationary_teeth: f64) -> Vec<f64> {
    let offsets = parse_values(arg_matches, name);
    for offset in offsets.iter() {
        if *offset < 0.0 {
            panic!("{}: {} is negative", name, offset);
        }
        if *offset >= stationary_teeth {
            panic!(
                "{}: {} is not less than stationary-teeth: {}",
                name, offset, stationary_teeth
            );
        }
    }
    offsets
}

fn parse_scale(arg_matches: &ArgMatches, name: &str) -> Scale {
    match arg_matches.value_of(name).unwrap() {
        "same" => Scale::Same,
        value if value.starts_with('x') => Scale::Factor(
            value[1..]
                .parse::<f64>()
                .unwrap_or_else(|error| panic!("{} in {}: {}", error, name, value)),
        ),
        value if value.ends_with("pt") => Scale::Size(
            value[..(value.len() - 2)]
                .parse::<f64>()
                .unwrap_or_else(|error| panic!("{} in {}: {}", error, name, value)),
        ),
        value => panic!("invalid {}: {}", name, value),
    }
}

fn parse_polyline(arg_matches: &ArgMatches, name: &str) -> Polyline {
    let value = arg_matches.value_of(name).unwrap();
    if let Result::Ok(sides) = value.parse::<usize>() {
        if sides == 0 {
            panic!("{} sides are zero", name);
        }
        regular_polyline(sides)
    } else {
        load_polyline_from_svg_file(value)
    }
}

fn regular_polyline(sides: usize) -> Polyline {
    if sides == 1 {
        parse_svg(
            "\
            <svg xmlns='http://www.w3.org/2000/svg'>\
            <path d='\
            M 1, 1\
            m -1, 0\
            a 1,1 0 1,0 2,0\
            a 1,1 0 1,0 -2,0\
            '/>\
            </svg>\
            ",
        )
        .unwrap()
        .pop()
        .unwrap()
    } else {
        let mut polyline = vec![Point { x: 0.0, y: 0.0 }; sides];
        let angle = (2 as f64) * PI / (sides as f64);
        for side in 0..sides {
            polyline[side].x = (angle * side as f64).sin();
            polyline[side].y = (angle * side as f64).cos();
        }
        polyline
    }
}
