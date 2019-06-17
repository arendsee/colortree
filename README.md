# colortree

## Installation

Installing from github requires the following command in the R shell:

```R
devtools::install_github("arendsee/colortree")
```


Or you can install locally if you want to edit the code

```bash
git clone https://github.com/arendsee/colortree
cd colortree
R
```

In the R shell

```R
R> devtools::document()
R> devtools::install()
```

## Using the script

The executable file is the `inst/colortree.R`. So you can symlink this to your
path. For example:

```bash
ln -s $PWD/inst/colortree.R $HOME/bin/colortree
```

Then you can edit the colortree.R file directly to customize your scripts.


## Execution

If the executable is linked to PATH as above, you can call it as so:

```bash
colortree mytree.tre
```

Which will product the file `mytree.pdf`.
