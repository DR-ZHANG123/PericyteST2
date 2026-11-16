# 可复现性工具（R 版）：读配置、固定种子、写 results/<stage>/PROVENANCE.json。
suppressPackageStartupMessages({library(yaml); library(jsonlite); library(tools)})

p22_root <- function() {
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "config", "run.yaml"))) {
    up <- dirname(d); if (up == d) stop("找不到 config/run.yaml，请在项目目录下运行"); d <- up
  }
  d
}

load_config <- function() yaml::read_yaml(file.path(p22_root(), "config", "run.yaml"))

set_global_seed <- function(seed) { set.seed(seed); invisible(seed) }

.relpath <- function(p, root) sub(paste0("^", root, "/"), "", normalizePath(p, mustWork = FALSE))

write_provenance <- function(stage, inputs, outputs, seed = NULL, extra = list()) {
  root <- p22_root()
  out_dir <- file.path(root, "results", stage); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cfg <- file.path(root, "config", "run.yaml")
  md5s <- function(v) { v <- v[file.exists(v)]; setNames(as.list(unname(md5sum(v))), vapply(v, .relpath, "", root = root)) }
  commit <- tryCatch(system2("git", c("-C", root, "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
                     error = function(e) "unknown", warning = function(w) "unknown")
  pk <- c("limma", "edgeR", "sva", "fgsea", "msigdbr")
  pk <- pk[vapply(pk, requireNamespace, logical(1), quietly = TRUE)]
  rec <- c(list(stage = stage,
                finished_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
                seed = seed,
                config_sha = substr(digest_sha256(cfg), 1, 16),
                git_commit = if (length(commit)) commit[1] else "unknown",
                r_version = R.version.string,
                conda_env = basename(dirname(dirname(R.home()))),   # R.home() = <env>/lib/R
                packages = setNames(as.list(vapply(pk, function(x) as.character(packageVersion(x)), "")), pk),
                inputs = md5s(inputs), outputs = md5s(outputs)), extra)
  path <- file.path(out_dir, "PROVENANCE.json")
  writeLines(toJSON(rec, auto_unbox = TRUE, pretty = TRUE, null = "null"), path)
  invisible(path)
}

digest_sha256 <- function(f) {
  out <- system2("sha256sum", f, stdout = TRUE)
  strsplit(out, " ")[[1]][1]
}
