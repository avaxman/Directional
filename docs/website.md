# Building the Website

For developers who want to contribute to the website/documentation of Directional.
If you want to preview changes to the Directional website before a commit, you can follow the instructions below.

### Prerequisites

1. If you do not already have it, install conda on your machine. We recommend using [miniconda3](https://docs.conda.io/en/latest/miniconda.html). On Linux you can run:
    ```bash
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh
    ```
2. Install the conda environment for the website:
   ```
   conda env create -f directional-website.yml
   ```

### Using mkdocs

1. Activate the conda environment installed on the previous step:
   ```bash
   conda activate directional-website
   ```
2. Preview the website locally
   (run this command in the root folder of the Directional repository):
   ```bash
   mkdocs serve
   ```

!!! tip

Dead links can be checked using the [LinkChecker](https://linkchecker.github.io/linkchecker/)
tool. Run the website locally, then run LinkChecker on it:
```bash
linkchecker http://127.0.0.1:8000
```
!!! note
    The reason we are using `python -m mkdocs serve` instead of `mkdocs serve` directly is because we are using local extensions for mkdocs. Those extensions are located in the `scripts/` folder of libigl. Running `mkdocs` as a module adds the current directory to the `PYTHONPATH`, allowing us to load those extensions without installing them on the system or in a virtualenv.

## References

- [MkDocs](http://www.mkdocs.org/)
- [Material Theme](https://squidfunk.github.io/mkdocs-material/)
