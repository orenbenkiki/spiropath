//! Functions for implementing the executable spiropath program.

use crate::library::load_polyline_from_svg_file;
use crate::library::print_svg_polylines;
use crate::library::spiropath;
use crate::library::Location;
use crate::library::Options;
use crate::library::Scale;
use clap::crate_version;
use clap::App;
use clap::Arg;
use clap::ArgMatches;
use std::f64::consts::PI;
use std::fs::File;
use std::io::stdout;
use std::io::BufWriter;
use svg2polylines::CoordinatePair as Point;
use svg2polylines::Polyline;

/// A complete main function for the spiropath program.
pub fn main(flags: &[String]) {
    let arg_matches = app().get_matches_from(flags);

    let mut options = arg_options(&arg_matches);
    options.validate();

    let stationary = parse_polyline(
        &arg_matches,
        "stationary",
        options.stationary_teeth,
        options.tolerance,
    );

    let rotating = parse_polyline(
        &arg_matches,
        "rotating",
        options.rotating_teeth,
        options.tolerance,
    );

    let polylines = spiropath(stationary, rotating, options);

    let output_path = arg_matches.value_of("output").unwrap();
    if output_path == "-" {
        print_svg_polylines(&polylines, &mut BufWriter::new(stdout())); // NOT TESTED
    } else {
        print_svg_polylines(
            &polylines,
            &mut BufWriter::new(
                File::create(output_path)
                    .unwrap_or_else(|error| panic!("{} creating output: {}", error, output_path)),
            ),
        );
    };
}

fn app() -> App<'static, 'static> {
    App::new("spiropath")
        .about("\nGeneralized Spirograph using arbitrary SVG paths.")
        .after_help(
            "\
            PROCESS:\n\
            \n\
            - Take as input a stationary SVG path and a rotating SVG path;\n  \
              due to svg2polylines limitations, arcs are not supported.\n  \
              Or, generate a circle or a regular polygon instead.\n\
            \n\
            - Scale these so their total length is an integer,\n  \
              equivalent to a number of virtual gear \"teeth\".\n\
            \n\
            - Pick the top-most point of the stationary path,\n  \
              and the top-most point of the (rotated) rotating path.\n  \
              Position the rotating path so its top point is at some offset(s)\n  \
              relative to the top-point of the stationary path.\n\
            \n\
            - Rotate the rotating path around the stationary path,\n  \
              until it returns to its start position.\n  \
              Ignore collisions; consider only the contact point.\n\
            \n\
            - Generate the path traced by a rotating point,\n  \
              which is fixed relative to the rotating path.\n  \
              Allow this point to be outside the rotating path.\n\
            \n\
            - Optionally also include the stationary path.\n\
            \n\
            - Scale the result; coordinates are in SVG \"pt\" units.\n\
            \n\
            - Print this as SVG file.\
            ",
        )
        .version(crate_version!())
        .arg(
            Arg::with_name("stationary")
                .long("stationary")
                .short("s")
                .value_name("FILE or CIRCLE or COUNT")
                .help(
                    "SVG file containing the rotating path,\n\
                       CIRCLE for generating a radius-100 circle,\n\
                       or the number of sides for a radius-100 polygon:\n\
                       2 - line, 3 - triangle, 4 - square, etc.\n",
                )
                .default_value("CIRCLE"),
        )
        .arg(
            Arg::with_name("rotating")
                .long("rotating")
                .short("r")
                .value_name("FILE or CIRCLE or COUNT")
                .help(
                    "SVG file containing the rotating path,\n\
                       CIRCLE for generating a radius-100 circle,\n\
                       or the number of sides for a radius-100 polygon:\n\
                       2 - line, 3 - triangle, 4 - square, etc.\n",
                )
                .default_value("CIRCLE"),
        )
        .arg(
            Arg::with_name("output")
                .long("output")
                .short("o")
                .value_name("FILE")
                .help("SVG file to write the output spiropath into;\nspecify \"-\" for STDOUT")
                .default_value("-"),
        )
        .arg(
            Arg::with_name("location")
                .long("location")
                .short("L")
                .value_name("SIDE")
                .possible_values(&["inside", "outside"])
                .help("Location of the rotating path relative to the stationary path\n")
                .default_value("outside"),
        )
        .arg(
            Arg::with_name("stationary-teeth")
                .long("stationary-teeth")
                .short("S")
                .value_name("COUNT")
                .help(
                    "Number of virtual teeth on (total length of)\n\
                      the stationary path",
                )
                .default_value("47"),
        )
        .arg(
            Arg::with_name("rotating-teeth")
                .long("rotating-teeth")
                .short("R")
                .value_name("COUNT")
                .help(
                    "Number of virtual teeth on (total length of)\n\
                      the rotating path",
                )
                .default_value("7"),
        )
        .arg(
            Arg::with_name("tolerance")
                .long("tolerance")
                .short("T")
                .value_name("FRACTION")
                .help(
                    "Linear path approximation tolerance,\n\
                      as a fraction of the teeth size",
                )
                .default_value("0.001"),
        )
        .arg(
            Arg::with_name("point")
                .long("point")
                .short("P")
                .help(
                    "Coordinates of the traced point,\n\
                       relative to the raw rotating path\n",
                )
                .value_name("X,Y")
                .number_of_values(2)
                .require_delimiter(true)
                .default_value("100,0"),
        )
        .arg(
            Arg::with_name("angle")
                .long("angle")
                .short("A")
                .value_name("DEGREES")
                .help(
                    "Angle to use for picking the top point of the rotating path,\n\
                      in degrees",
                )
                .default_value("0.0"),
        )
        .arg(
            Arg::with_name("offset")
                .long("offset")
                .short("O")
                .value_name("FRACTION(S)")
                .help(
                    "Offset(s) of start position(s) of rotating path,\n\
                      relative to top point of stationary path,\n\
                      as a fraction of the teeth size;\n\
                      repeat for including multiple paths in the output\n",
                )
                .multiple(true)
                .use_delimiter(true)
                .min_values(1)
                .default_value("0.0"),
        )
        .arg(
            Arg::with_name("x-scale")
                .long("x-scale")
                .short("X")
                .value_name("SCALE")
                .help(
                    "Scaling of the output SVG, one of:\n\
                       \"<size>pt\" / \"x<factor>\" / \"same\" as y-scale\n",
                )
                .default_value("x1.0"),
        )
        .arg(
            Arg::with_name("y-scale")
                .long("y-scale")
                .short("Y")
                .value_name("SCALE")
                .help(
                    "Scaling of the output SVG, one of:\n\
                       \"<size>pt\" / \"x<factor>\" / \"same\" as x-scale\n",
                )
                .default_value("x1.0"),
        )
        .arg(
            Arg::with_name("mirror-rotating")
                .long("mirror-rotating")
                .short("M")
                .help("Mirror the rotating path"),
        )
        .arg(
            Arg::with_name("include-stationary")
                .long("include-stationary")
                .short("I")
                .help("Include the stationary path in the output"),
        )
}

fn arg_options(arg_matches: &ArgMatches) -> Options {
    let location = parse_location(arg_matches, "location");
    let stationary_teeth = parse_count(arg_matches, "stationary-teeth");
    let rotating_teeth = parse_count(arg_matches, "rotating-teeth");
    let tolerance = parse_fraction(arg_matches, "tolerance");
    let traced_point = parse_point(arg_matches, "point");
    let initial_rotating_angle = parse_value(arg_matches, "angle");
    let mirror_rotating = arg_matches.is_present("mirror-rotating");
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
        mirror_rotating,
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
            // BEGIN NOT TESTED
            panic!(
                "{} in {}: {}",
                error,
                name,
                arg_matches.value_of(name).unwrap()
            )
            // END NOT TESTED
        });
    if value == 0 {
        panic!("{} is zero", name); // NOT TESTED
    }
    value
}

fn parse_fraction(arg_matches: &ArgMatches, name: &str) -> f64 {
    let value = parse_value(arg_matches, name);
    if value <= 0.0 {
        panic!("{}: {} is not positive", name, value); // NOT TESTED
    }
    value
}

fn parse_value(arg_matches: &ArgMatches, name: &str) -> f64 {
    arg_matches
        .value_of(name)
        .unwrap()
        .parse::<f64>()
        .unwrap_or_else(|error| {
            // BEGIN NOT TESTED
            panic!(
                "{} in {}: {}",
                error,
                name,
                arg_matches.value_of(name).unwrap()
            )
            // END NOT TESTED
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
    let mut offsets = parse_values(arg_matches, name);
    for offset in offsets.iter_mut() {
        while *offset < 0.0 {
            *offset += stationary_teeth; // NOT TESTED
        }
        while *offset > stationary_teeth {
            *offset -= stationary_teeth; // NOT TESTED
        }
        if *offset >= stationary_teeth {
            // BEGIN NOT TESTED
            panic!(
                "{}: {} is not less than stationary-teeth: {}",
                name, offset, stationary_teeth
            );
            // END NOT TESTED
        }
    }
    offsets
}

fn parse_scale(arg_matches: &ArgMatches, name: &str) -> Scale {
    match arg_matches.value_of(name).unwrap() {
        "same" => Scale::Same,
        value if value.starts_with('x') => Scale::Factor(
            // BEGIN NOT TESTED
            value[1..]
                .parse::<f64>()
                .unwrap_or_else(|error| panic!("{} in {}: {}", error, name, value)),
            // END NOT TESTED
        ),
        value if value.ends_with("pt") => Scale::Size(
            value[..(value.len() - 2)]
                .parse::<f64>()
                .unwrap_or_else(|error| panic!("{} in {}: {}", error, name, value)),
        ),
        value => panic!("invalid {}: {}", name, value), // NOT TESTED
    }
}

fn parse_polyline(arg_matches: &ArgMatches, name: &str, teeth: usize, tolerance: f64) -> Polyline {
    let value = arg_matches.value_of(name).unwrap();
    if value == "CIRCLE" {
        // BEGIN NOT TESTED
        let max_step_angle = 2.0 * (1.0 - tolerance / 100.0).acos();
        let minimal_teeth = 2.0 * PI / max_step_angle;
        let factor = (minimal_teeth / teeth as f64).ceil();
        regular_polyline(teeth * factor as usize)
        // END NOT TESTED
    } else if let Result::Ok(sides) = value.parse::<usize>() {
        if sides < 2 {
            panic!("{} sides: {} are less than 2", name, sides); // NOT TESTED
        }
        regular_polyline(sides)
    } else {
        load_polyline_from_svg_file(value, tolerance) // NOT TESTED
    }
}

fn regular_polyline(sides: usize) -> Polyline {
    let angle = 2.0 * PI / (sides as f64);
    let mut polyline = vec![Point { x: 0.0, y: 0.0 }; sides];
    for (side, point) in polyline.iter_mut().enumerate() {
        point.x = 100.0 * (angle * side as f64).sin();
        point.y = 100.0 * (angle * side as f64).cos();
    }
    polyline
}
