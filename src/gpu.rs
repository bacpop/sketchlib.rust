
#[link(name = "vector_add", kind = "static")]
extern "C" {
    fn vectorAdd_main();
}

fn call_gpu() {
    unsafe {
        vectorAdd_main();
    }
}