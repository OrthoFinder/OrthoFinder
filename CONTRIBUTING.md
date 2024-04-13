# Contributing

Contributions are welcome, and they are greatly appreciated! Every little bit
helps, and credit will always be given.

## Types of Contributions

### Report Bugs

If you are reporting a bug, please post an [Issue](https://github.com/OrthoFinder/OrthoFinder/issues). Please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

### Fix Bugs

Look through the [OrthoFinder Issues log](https://github.com/OrthoFinder/OrthoFinder/issues) for bugs. Anything tagged with "bug" 
and "help wanted" is open to whoever wants to implement it.

### Implement Features

Look through the [OrthoFinder Issues log](https://github.com/OrthoFinder/OrthoFinder/issues) for features. Anything tagged with 
"enhancement" and "help wanted" is open to whoever wants to implement it.

### Write Documentation

You can never have enough documentation! Please feel free to contribute to any
part of the documentation, such as the official docs, docstrings, or even
on the web in blog posts, articles, and such.

### Submit Feedback

If you are proposing a feature:

* Set up an Issue or WIP Pull Request to kick-off discussions
* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a passion-driven project, and that contributions
  are welcome &#x1F604;

## Get Started!

Ready to contribute? Here's how to set up `system-identification` for local development.

1. Clone [OrthoFinder](https://github.com/OrthoFinder/OrthoFinder) locally
  
2. Install `OrthoFinder` in editable mode using `pip`:

    ```console
    $ python -m pip install -e .
    ```

3. Use `git` (or similar) to create a branch for local development and make your changes:

    ```console
    $ git checkout -b name-of-your-bugfix-or-feature
    ```

4. Push your branch and set up a new Pull Request: use `WIP` prefix to indicate 
   this is a work in progress (this way other collaborators have visibility of 
   your work even when it is unfinished)

5. Make commits and push them to the OrthoFinder GitHub repository when convenient. Please note each 
   time a batch of commits is pushed, OrthoFinder Git Actions will run the test suite

6. When you're done making changes, check that your changes conform to any code 
   formatting requirements that the tests pass. You should then add `OrthoFinder` as Reviewer, to obtain 
   QA-review and approval of your work

7. Once approval is provided, your branch will be merged

## Pull Request Guidelines

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include additional tests if appropriate.
2. If the pull request adds functionality, the docs should be updated.
3. The pull request should work for all currently supported operating systems and versions of Python.

## Code of Conduct

Please note that the `OrthoFinder` project is released with a Code of Conduct. By contributing to this project you agree to abide by its terms.

## Commiting Guidelines

`OrthoFinder` uses [Python Semantic Release (PSR)](https://python-semantic-release.readthedocs.io/en/latest/) to automatically bump version numbers based on keywords it finds in commit messages.

The idea of PSR is to use a standardized commit message format and syntax, which PSR can parse to determine how to increment the version number. The default commit message format used by PSR is the [Angular commit style](https://github.com/angular/angular.js/blob/master/DEVELOPERS.md#commit-message-format), which looks like this:
```console
<type>(optional scope): short summary in present tense

(optional body: explains motivation for the change)

(optional footer: note BREAKING CHANGES here, and issues to be closed)
```

`<type>` refers to the kind of change made and is usually one of:

+ **feat**: A new feature. 
+ **fix**: A bug fix.
+ **docs**: Documentation changes.
+ **style**: Changes that do not affect the meaning of the code (white-space, formatting, missing semi-colons, etc).
+ **refactor**: A code change that neither fixes a bug nor adds a feature.
+ **perf**: A code change that improves performance.
+ **test**: Changes to the test framework.
+ **build**: Changes to the build process or tools.

`scope` is an optional keyword that provides context for where the change was made. It can be anything relevant to your package or development workflow (e.g., it could be the module or function name affected by the change).

Different text in the commit message will trigger PSR to make different kinds of releases:

+ A `<type>` of fix triggers a patch version bump, e.g.:
```console
$ git commit -m "fix(notes): fix confusing message in notes"
```
+ A `<type>` of feat triggers a minor version bump, e.g.:

```console
$ git commit -m "feat(package): add scoring matrix option to the package"
```
+ The text `BREAKING CHANGE:` in the footer will trigger a major release, e.g.:
```console
$ git commit -m "build(package): reconfig the package using src layout
$ 
$ BREAKING CHANGE: the scripts_of folder does not exit anymore."
```


