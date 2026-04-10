use std::path::Path;

use nerpa_ms_core::testing::dp_vs_bruteforce::test_dp_vs_bruteforce;

#[test]
fn dp_vs_bruteforce_curated() {
    let tests_json = Path::new("tests/fixtures/dp_vs_bruteforce_tests_curated.json");
    let out_dir = Path::new("test_output/dp_vs_bruteforce_svgs"); 
    test_dp_vs_bruteforce(tests_json, out_dir);
}

//#[test]
fn dp_vs_bruteforce_rand() {
    let tests_json = Path::new("tests/fixtures/dp_vs_bruteforce_tests_rand.json");
    let out_dir = Path::new("test_output/dp_vs_bruteforce_svgs"); 
    test_dp_vs_bruteforce(tests_json, out_dir);
}

