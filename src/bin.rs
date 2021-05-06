// FILE MAYBE TESTED

use std::env::args;

fn main() {
    let flags: Vec<String> = args().collect();
    spiropath::program::main(&flags);
}
