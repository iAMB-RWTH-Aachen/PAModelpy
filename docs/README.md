# Website

This website is built using [Docusaurus](https://docusaurus.io/), a modern static website generator.
You will need NodeJS to be able to run Docusaurus for development.

### Installation
To use Docusaurus, install NodeJS. You will find a comprehensive guide [here](https://nodejs.org/en/download/package-manager#debian-and-ubuntu-based-linux-distributions) for linux systems.

```
$ npm install
```

### Local Development

```
$ npm run start
```

This command starts a local development server and opens up a browser window. Most changes are reflected live without having to restart the server.

### Build

```
$ npm run build
```

This command generates static content into the `build` directory and can be served using any static contents hosting service.

### Deployment

Using SSH:

```
$ USE_SSH=true npm run deploy
```

Not using SSH:

```
$ GIT_USER=<Your GitHub username> npm run deploy
```

If you are using GitHub pages for hosting, this command is a convenient way to build the website and push to the `gh-pages` branch.

# Documenting your code

## Setup

This installation makes use of the Python module "pydoc-markdown". If you want to use it, install this module and execute it in the base directory of this repo. The provided pydoc-markdown.yml file contains all necessary configurations needed. It will generate the API reference from the docstrings in the PAModelpy package. Please note, that you have to use Google docstring format if you want to make changes to the PAModelpy project.
