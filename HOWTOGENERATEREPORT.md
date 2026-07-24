# Lipidomics

## Générer le rapport

``` bash
cd /Users/loictalignani/research/project/lipidomics/analysis \&\& Rscript -e "rmarkdown::render('Lipidomics\_Analysis\_Report.Rmd', output\_format = 'html\_document')

Rscript -e "rmarkdown::render('Lipidomics\_Analysis\_Report\_D3.Rmd', output\_format = 'html\_document')"
```

# Si vous voulez le générer depuis R directement

rmarkdown::render("Lipidomics_Analysis_Report.Rmd")

# Pour un format PDF (nécessite LaTeX)

rmarkdown::render("Lipidomics_Analysis_Report.Rmd", output_format = "pdf_document")

# Pour un Word document

rmarkdown::render("Lipidomics_Analysis_Report.Rmd", output_format = "word_document")
