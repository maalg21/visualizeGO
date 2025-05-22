# visualizeGO

**visualizeGO** is an R package for visualising, analysing, and summarising gene ontology (GO) term enrichment results. It allows graphical representation of GO term clusters based on semantic similarity or co-occurrence networks, facilitating their interpretation.

---

## 🚀 Instalation

You can install the development version from GitHub with:

```r
# Install 'devtools' if you don't already have it
install.packages("devtools")

# Instala visualizeGO desde GitHub
devtools::install_github("tu_usuario/visualizeGO")
```
## ⚡ Quick use

``` r
library(visualizeGO)

# Assuming you have an object ‘go_results’ (from enrichGO, gProfiler, etc.)
clustered <- clusterGO(go_results)
visualize_clusters(clustered)

```

## 📘 Complete tutorial

See the vignette for a complete example of the workflow, including preprocessing, clustering, visualisation, and exporting results:
``` r
vignette("visualizeGO-intro")
```

## 📬 Contact

If you have any questions, suggestions, or would like to contribute to the development of **visualizeGO**, please do not hesitate to contact me:

- 📧 Email: [maalg@unileon.es](mailto:maalg@unileon.es)
- 🐙 GitHub Issues: [Open an issue](https://github.com/usuario/visualizeGO/issues)
- 🧪 Report bugs or request improvements: use the "Issues"" tab in the repository.

Thanks for using `visualizeGO`!
