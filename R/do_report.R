columns <- list()

columns$bootstrap <- c("Gene Size",
                       "Reference Connectivity",
                       "Test Connectivity",
                       "MDC Unadusted",
                       "MDC Adjusted",
                       "MDC Stat",
                       "MDC Stdev",
                       "MDC P-Val",
                       "MDC FDR",
                       "KS Stat",
                       "KS Stdev",
                       "KS P-Val",
                       "KS FDR")

columns$permutation <- c("Gene Size",
                         "Reference Connectivity",
                         "Test Connectivity",
                         "MDC Unadusted",
                         "MDC Adjusted",
                         "MDC P-Val",
                         "MDC FDR",
                         "KS P-Val",
                         "KS FDR")

#' @title Temporary
#' @description 
#' @param placeholder A placeholder
#' @return A return
#'
#' @importFrom magrittr %>%
#'
#' @export
do_report <- function(output, directory) {
    require(knitr)

    # Create main directory if not exists
    if (!file.exists(directory)) {
        dir.create(directory)
    }

    # Create folder for PDF plots if not exists
    pdf_dir <- file.path(directory, "PDFs")
    if (!file.exists(pdf_dir)) {
        dir.create(pdf_dir)
    }

    # Create folder for CDF plots if not exists
    cdf_dir <- file.path(directory, "CDFs")
    if (!file.exists(cdf_dir)) {
        dir.create(cdf_dir)
    }

    # Create Report File
    report_file <- file.path(directory, "report.rmd")
    head <- c("---", "title: 'Differential Connectivity'", "---")
    style <- c("<style type='text/css'>", "table {", "font-size: 8px;", "}", "</style>")

    mdc <- as.data.frame(output$info$mod_n)

    mdc <- cbind(mdc, as.data.frame(output$stat[c("mods_mc_r",
                                                  "mods_mc_t",
                                                  "mods_mdc_org",
                                                  "mods_mdc_adj")]))
    mdc <- cbind(mdc, output$iter)

    if (ncol(mdc) == 9) {
        colnames(mdc) <- columns$permutation
    } else {
        colnames(mdc) <- columns$bootstrap
    }
    mod_names <- rownames(mdc)
    write.table(mdc, file.path(directory, "mdc.txt"),
                quote=F,
                col.names=T,
                row.names=T,
                sep="\t")

    mdc_link <- paste("[Tab-Delimited Table](", "mdc.txt", ")", sep = "")

    if (!is.null(output$bg)){
        bg <- as.data.frame(output$bg[c("mc_r_bg", "mc_t_bg")])
        colnames(bg) <- c("Reference Connectivity", "Test Connectivity")
        write.table(bg, file.path(directory, "bg.txt"),
                    quote=F,
                    col.names=T,
                    row.names=F,
                    sep="\t")
        bg_link <- paste("[Tab-Delimited Table](", "bg.txt", ")", sep = "")
    }

    cat("Saving PDF/CDF plots for each module...\n")

    # Create PDF files
    lapply(mod_names, function(x) {
        png(file.path(pdf_dir, paste(x, "_pdf.png", sep = "")), height = 400, width = 700)
        print(output$plots$plot_pdf[[x]])
        dev.off()
    })

    # Create CDF files
    lapply(mod_names, function(x) {
        png(file.path(cdf_dir, paste(x, "_cdf.png", sep = "")), height = 400, width = 700)
        print(output$plots$plot_cdf[[x]])
        dev.off()
    })

    mdc$pdf <- paste("[![](", paste(file.path("PDFs", paste(mod_names, "_pdf.png", sep = "")), sep = ""), ")](",paste(file.path("PDFs", paste(mod_names, "_pdf.png", sep = "")),")", sep = ""))
    mdc$cdf <- paste("[![](", paste(file.path("CDFs", paste(mod_names, "_cdf.png", sep = "")), sep = ""), ")](",paste(file.path("CDFs", paste(mod_names, "_cdf.png", sep = "")),")", sep = ""))

    # Create markdown table
    mdc_rmd <- kable(mdc, digits = 4, align = "l")

    # Create report
    cat("Creating report ... \n")

    if (exists("bg")) {
        bg_rmd <- kable(bg, digits = 4, align = "l")
        writeLines(c(head, style, "### Background Statistics", bg_link, " ", bg_rmd, "### MDC Statistics", mdc_link, " ", mdc_rmd), report_file)
    } else {
        writeLines(c(head,style, "### MDC Statistics", mdc_link, " ", mdc_rmd), report_file)
    }
    rmarkdown::render(report_file, quiet=T)
    cat("Success: Report is saved as ", file.path(directory, "report.html"))
}
