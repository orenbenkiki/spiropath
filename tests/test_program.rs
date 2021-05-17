use spiropath::program::main;
use std::fs::create_dir_all;
use std::fs::read;
use std::fs::write;

#[macro_export]
macro_rules! test_name {
    () => {{
        fn f() {}
        fn type_name_of<T>(_: T) -> &'static str {
            std::any::type_name::<T>()
        }
        let name = type_name_of(f);
        let prefix = &name[..name.len() - 3];
        let offset = prefix.rfind("::").unwrap();
        &prefix[offset + 2..]
    }};
}

#[macro_export]
macro_rules! test_case {
    ($name:ident, $suffix:literal, $flags:expr) => {
        #[test]
        fn $name() {
            for directory in vec!["tests/expected", "tests/actual"].iter() {
                create_dir_all(directory).unwrap_or_else(|_| {
                    // BEGIN NOT TESTED
                    panic!("failed to create {} results directory", directory)
                    // END NOT TESTED
                });
            }

            let file_name = format!("{}.{}", test_name!(), $suffix);
            let mut flags: Vec<String> = $flags.iter().map(|string| string.to_string()).collect();
            flags.push("--output".to_string());
            flags.push(format!("tests/actual/{}", file_name));
            main(&flags);
            impl_assert_output(&file_name);
        }
    };
}

fn impl_assert_output(file_name: &str) {
    let actual_path = format!("tests/actual/{}", file_name);
    let actual_bytes = read(actual_path.clone()).unwrap();

    let expected_path = format!("tests/expected/{}", file_name);
    let expected_bytes = read(expected_path.clone()).unwrap_or_else(|_| {
        // BEGIN NOT TESTED
        write(expected_path.clone(), &actual_bytes).unwrap_or_else(|_| {
            panic!("failed to write expected results file {}", expected_path);
        });
        eprintln!(
            "WARNING: created expected results file {}, verify its contents",
            expected_path
        );
        actual_bytes.clone().to_vec()
        // END NOT TESTED
    });

    assert!(
        expected_bytes == actual_bytes,
        "The actual results file {} is different from the expected results file {}",
        expected_path,
        actual_path
    );
}

test_case! {
    square,
    "svg",
    vec!["test", "-r", "4", "-s", "4", "-R", "2", "-S", "8", "-I", "-X", "100pt", "-Y", "100pt"]

}
