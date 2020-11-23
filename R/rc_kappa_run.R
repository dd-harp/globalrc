if (dir.exists("gen_scaled_ar")) {
    source(file.path("gen_scaled_ar", "rc_kappa.R"))
} else {
    source(file.path(".", "rc_kappa.R"))
}

main()
