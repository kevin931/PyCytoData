# PythonTemplate
> A template repository for a Python package

If you're looking for a quick template for your python project, you're at the right palce! This
is a skeleton project to kickstart (or Press Button Start) your next Python Project on the
horizon! Simply use the project without worrying about anything else!

For documentation and unittesting, you don't have to use the same framework as I'm using. If
you don't like it, leave parts of the folder out! No worry: I won't be offendeds!

## Packages You Will Need

You will need some packages to get started! I recommend using a ``conda`` environment, but
you can use whatever virtual environment framework you want!

### Python

Duh! You will need a relatively modern (I'm talking 3.7 and up) althought the ``setup.py`` specifies
3.9.

### For Documentation

I use ``sphinx`` and ``readthedocs`` for documentation! You will need to install the following packages:

```
sphinx
sphinx-rtd-theme
sphinx-git
sphinxcontrib-autoprogram
sphinx-autodoc-typehints
```

### For Unit Testing

You will need the following packages:

```
pytest
pytest-cov
pytest-mock
coverage
```

Technically, only ``pytest`` is strictly necessary if you don't want to use coverage or ``mocker``. 

## Documentation

I use ``sphinx`` for documentation, which is super easy (well, it takes a little bit of time)! A lot of the 
``rst`` files I added myself, but to get the basic structure, you will need to run:

```shell
cd docs
sphinx-quickstart
```
This will set up your directory with the proper files, just like mine, and you can copy some of the ``rst`` there!

### I'm a Newbie, what should I do?

I recommend you look through [this guide](https://docs.readthedocs.io/en/stable/intro/getting-started-with-sphinx.html)

### What does documentation look like?

Here is the link to how this [repository's documentation](pythontemplateproject.readthedocs.io)

### How do I build the documentation locally?

Easy! [Look here!](https://pythontemplateproject.readthedocs.io/en/latest/change/build.html)


## Unit Tests

[Here is a tutorial](https://pythontemplateproject.readthedocs.io/en/latest/change/build.html) on how to run tests! You may
also want to look into ``pytest`` for more details!

## What do I need to change?

You will need to change a number of things! Don't panic. They're easy!

### Directory and Python

- Change the ``PythonTemplate`` directory to your project name!
- Put in your own package's details in ``setup.py``.

### External Services

- Create a public GitHub Repo
- Set up documentation at [ReadTheDocs](https://readthedocs.org/) by linking your repository
- Register and link your GitHub account at CodeCov for coverage statistics

### Documentation
- Set up the documentation with ``sphinx`` as documented above. 
- Go through the documentation in the ``docs`` directory and add and change accordingly.

### GitHub Actions
- Change ``line 27`` in ``.github/workflows/ci.yaml`` to your own package name

