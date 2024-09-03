/*
   This installs R packages that are not included in the
   base tidyverse singularity image we are using.

   This way of doing this decreases reproducibility of 
   the pipeline.  It uses install.package() to install
   whatever the latest version of the corresponding
   packages. So not a fixed specific version.   
   But it's not straightforward to install
   specific versions of R packages using install.package()

   An alternative approach would be to make or find a mulled 
   (multi-package) singularity image with tidyverse and
   the packages we want.

   You could also try installing from tarballs, but you run
   into dependency rabbit-hole issues.

   See:
   https://stackoverflow.com/questions/17082341/installing-older-version-of-r-package
   https://lpembleton.rbind.io/posts/mulled-biocontainers/
   https://github.com/rocker-org/rocker-versioned2/blob/master/dockerfiles/tidyverse_devel.Dockerfile
*/
process SETUP_R_DEPENDENCIES {
  label      'process_low'
  tag        "${R_packages}"
  publishDir "${params.outdir}"

  // singularity info for this process
  if (workflow.containerEngine == 'singularity'){
      container "docker://rocker/tidyverse:4.3.2"
  }

  input:
  val (R_packages)        // a space-separated string of R package names to install (e.g.: "rstatix ggpubr") 

  output:
  val  (true)           , emit: setup_complete // a flag to indicate this is complete
  path ("R_lib_dir")    , emit: R_lib_dir
  // path "versions.yml"   , emit: versions

  script:

  if (workflow.containerEngine == 'singularity') {
  """
    #!/usr/bin/env Rscript

    # split up package names into a vector
    packages = strsplit("${R_packages}", " ")[[1]]

    # will install packages into a new, local lib dir
    r_lib_dir = "R_lib_dir"

    # create the new R_lib_dir
    dir.create(r_lib_dir, showWarnings = FALSE)

    # install each package
    for (package in packages) {
      # install package into specified library directory
      install.packages(package, lib=r_lib_dir)
    
      # confirm it worked
      library(package, character.only=T, lib.loc=r_lib_dir)
    }

  """
  } else {
  """
    echo "setup not necessary for conda environment"
    # just create an empty directory 
    mkdir R_lib_dir

    touch versions.yml
  """
  }
}
