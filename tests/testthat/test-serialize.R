test_that("serializion can read what it wrote", {
  chunk1 = array(
    1:12, dim = c(2, 3, 2),
    dimnames = list(row=character(0), col=character(0), year=character(0))
    )
  attr(chunk1, "tile") <- c(2, 7)
  chunk2 = array(
    1:12, dim = c(2, 6),
    dimnames = list(row=3:4, col=6:11)
  )
  attr(chunk1, "tile") <- c(2, 8)
  test_chunks <- list(
    a = chunk1,
    b = chunk2
  )
  lapply(c("row", "col"), function(x) NULL)
  fn <- tempfile(fileext = ".h5")
  save_chunks(fn, test_chunks)
  cc <- load_chunks(fn)
  is.null(dimnames(cc$a)$row)
  is.null(dimnames(cc$a)$col)
  expect_true(all(dimnames(cc$b)$row == as.character(3:4)))
  unlink(fn)
})
