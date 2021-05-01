# spiropath
[![license](https://img.shields.io/crates/l/spiropath.svg)](https://github.com/orenbenkiki/spiropath/blob/master/LICENSE.txt)
[![crates.io](https://img.shields.io/crates/v/spiropath.svg)](https://crates.io/crates/spiropath)
[![Build Status](https://api.travis-ci.com/orenbenkiki/spiropath.svg?branch=master)](https://travis-ci.com/orenbenkiki/spiropath)
[![codecov](https://codecov.io/gh/orenbenkiki/spiropath/branch/master/graph/badge.svg)](https://codecov.io/gh/orenbenkiki/spiropath)
[![Docs](https://docs.rs/spiropath/badge.svg)](https://docs.rs/crate/spiropath)

Investigate the total state space of communicating finite state machines. Specifically, given a
model of a system comprising of multiple agents, where each agent is a non-deterministic state
machine, which responds to either time or receiving a message with one of some possible actions,
where each such action can change the agent state and/or send messages to other agents; Then the
code here will generate the total possible configurations space of the overall system, validate the
model for completeness, validate each system configuration for additional arbitrary correctness
criteria, and visualize the model in various ways.

## Installing

To install:

```
cargo install spiropath
```

## Using

This crate installs a `spiropath` program you can directly invoke from the command line: `spiropath
-s stationary.svg -r rotating.svg -o output.svg`, where `stationary.svg` contains a single SVG path
for the stationary shape and `rotating.svg` contains a single SVG path for the shape rotating around
it. The output SVG file will contain a single polyline tracing the path of a point affixed to the
rotating shape as it is moved around the rotating shape.

In addition this crate provides the same functionality (sans the file input/output) as a
`spiropath::spiropath` library function you can use in your own code.

A "classical" [Spirograph](https://en.wikipedia.org/wiki/Spirograph) uses only circles for both the
stationary and rotating shapes. Some Spirograph programs (or physical sets) provide more shapes but
in general these shapes are restricted to a combination of lines and arcs, and the overall
configuration (even in programs) is restricted by physics (e.g. the traced point must be internal to
the rotating shape).

In contrast, Spiropath places no restrictions on the stationary and rotating shapes (other than
requiring them to be continuous closed paths), or the position of the traced point: the rotating
shape is moved around the stationary one ignoring any collisions, and the traced point is moved with
it. You can therefore not only specify shapes using any cubic curves (rather than being restricted
to lines and arcs), but you can also use non-convex (or even self-intersecting) shapes, and in
general specify configurations which are impossible to achieve using a physical spirograph device
(e.g. placing the trace point outside the rotating shape).

In general the output SVG is expected to be further processed before the final visualization (e.g.
using [Inkscape](https://inkscape.org), so Spiropath doesn't provide options to control the path's
color, style, width, etc.

Spiropath does provide options to specify the relative sizes of the stationary and the rotating
shapes (the number of "teeth" along each of the shapes), the accuracy tolerance (as a fraction of
the tooth size), the size of the output (its bounding box), the initial rotation and offsets of the
rotating shape, the position of the traced point in the coordinates of the rotated shape, and
whether the rotating shape is placed "inside" or "outside" the stationary shape.

The resulting polyline will have a "lot" of points (depending on the tolerance); ideally we'd
convert it to a compact representation using SVG path elements such as arcs and curves. However, I
couldn't find a Rust library to do this. If you are aware of one, please let me know.

In the meanwhile, you can use external tools, for example using [Inkscape](https://inkscape.org/) to
either manually simplify the path as described in the "Simplification" section
[here](https://inkscape.org/doc/tutorials/advanced/tutorial-advanced.html) or automate this as
described
[here](https://stackoverflow.com/questions/7299564/automatizing-simplify-path-for-a-svg-file-using-inkscape).

## License

`spiropath` is distributed under the GNU AFFERO General Public License (Version 3.0). See the
[LICENSE](LICENSE.txt) for details.
