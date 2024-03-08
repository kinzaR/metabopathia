# Automate Package and Project Setup
# think about this : https://cran.r-project.org/web/packages/usethis/usethis.pdf
library(usethis)

# Create a new package -------------------------------------------------
path <- file.path(tempdir(), "metabopathia")
create_package(path)
#> ✔ Creating '/tmp/RtmpPqIkgo/mypkg/'
#> ✔ Setting active project to '/private/tmp/RtmpPqIkgo/mypkg'
#> ✔ Creating 'R/'
#> ✔ Writing 'DESCRIPTION'
#> Package: mypkg
#> Title: What the Package Does (One Line, Title Case)
#> Version: 0.0.0.9000
#> Authors@R (parsed):
#>     * First Last <first.last@example.com> [aut, cre] (YOUR-ORCID-ID)
#> Description: What the package does (one paragraph).
#> License: `use_mit_license()`, `use_gpl3_license()` or friends to pick a
#>     license
#> Encoding: UTF-8
#> Roxygen: list(markdown = TRUE)
#> RoxygenNote: 7.2.3
#> ✔ Writing 'NAMESPACE'
#> ✔ Setting active project to '<no active project>'
# only needed since this session isn't interactive
proj_activate(path)
#> ✔ Setting active project to '/private/tmp/RtmpPqIkgo/mypkg'
#> ✔ Changing working directory to '/tmp/RtmpPqIkgo/mypkg/'

# Modify the description ----------------------------------------------
use_mit_license("My Name")
#> ✔ Adding 'MIT + file LICENSE' to License
#> ✔ Writing 'LICENSE'
#> ✔ Writing 'LICENSE.md'
#> ✔ Adding '^LICENSE\\.md$' to '.Rbuildignore'

use_package("ggplot2", "Suggests")
#> ✔ Adding 'ggplot2' to Suggests field in DESCRIPTION
#> • Use `requireNamespace("ggplot2", quietly = TRUE)` to test if package is installed
#> • Then directly refer to functions with `ggplot2::fun()`

# Set up other files -------------------------------------------------
use_readme_md()
#> ✔ Writing 'README.md'
#> • Update 'README.md' to include installation instructions.

use_news_md()
#> ✔ Writing 'NEWS.md'

use_test("my-test")
#> ✔ Adding 'testthat' to Suggests field in DESCRIPTION
#> ✔ Adding '3' to Config/testthat/edition
#> ✔ Creating 'tests/testthat/'
#> ✔ Writing 'tests/testthat.R'
#> ✔ Writing 'tests/testthat/test-my-test.R'
#> • Edit 'tests/testthat/test-my-test.R'

x <- 1
y <- 2
use_data(x, y)
#> ✔ Adding 'R' to Depends field in DESCRIPTION
#> ✔ Creating 'data/'
#> ✔ Setting LazyData to 'true' in 'DESCRIPTION'
#> ✔ Saving 'x', 'y' to 'data/x.rda', 'data/y.rda'
#> • Document your data (see 'https://r-pkgs.org/data.html')

# Use git ------------------------------------------------------------
use_git()
#> ✔ Initialising Git repo
#> ✔ Adding '.Rproj.user', '.Rhistory', '.Rdata', '.httr-oauth', '.DS_Store', '.quarto' to '.gitignore'
