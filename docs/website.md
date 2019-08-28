# Building the Website

For developers who want to contribute to the website/documentation of Directional.
If you want to preview changes to the Directional website before a commit, you can follow the instructions below.

1. Install mkdocs and the material theme
   ```bash
   pip3 install -U --user mkdocs mkdocs-material
   ```
   Alternatively, use `pipenv` to install the dependencies:
   ```bash
   pip3 install -U --user pipenv
   pipenv install requests
   pipenv shell
   ```
2. Preview the website locally (in the root folder of the Directional project):
   ```bash
   python3 -m mkdocs serve
   ```
3. Build the website to generate the html locally (optional):
   ```bash
   python3 -m mkdocs build
   ```

!!! tip
    Dead links can be checked using the [LinkChecker](https://wummel.github.io/linkchecker/) tool. Run the website locally, then run LinkChecker on it:
    ```bash
    linkchecker http://127.0.0.1:8000
    ```

!!! note
    The reason we are using `python -m mkdocs serve` instead of `mkdocs serve` directly is because we are using local extensions for mkdocs. Those extensions are located in the `scripts/` folder of libigl. Running `mkdocs` as a module adds the current directory to the `PYTHONPATH`, allowing us to load those extensions without installing them on the system or in a virtualenv.

## TODO

- [ ] Use virtualenv to run a specific version of mkdocs/mkdocs-material, to avoid incompatible changes from breaking the build...

## References

- [MkDocs](http://www.mkdocs.org/)
- [Material Theme](https://squidfunk.github.io/mkdocs-material/)
