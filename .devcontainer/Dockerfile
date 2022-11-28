# references:
# https://github.com/revodavid/devcontainers-rstudio
# https://github.com/revodavid/devcontainers-rstudio/blob/main/.devcontainer/Dockerfile

# ------ Setup R 
# R version: 4, 4.1, 4.0, 4.2
ARG VARIANT="4.1"
#FROM rocker/tidyverse:${VARIANT}
# Includes packages: devtools tidyverse tidymodels glmnet randomForest xgboost 
FROM rocker/r-ver:${VARIANT}
# simpler, non full tidyverse install

# ------ Configure all additional options ------
# Use the [Option] comment to specify true/false arguments that should appear in VS Code UX
#
# [Option] Install zsh
ARG INSTALL_ZSH="false"
# [Option] Upgrade OS packages to their latest versions
ARG UPGRADE_PACKAGES="false"

RUN install2.r --error --skipinstalled --ncpus -1 \
    devtools \
    here \
    ggplot2 \
    remotes \
    data.table \
    readr \
    RCurl \
    && rm -rf /tmp/downloaded_packages

## Install summarytools and load some R package off the bat
#RUN R -e "devtools::install_github('https://github.com/dcomtois/summarytools/tree/0-8-9')"
#RUN R -e "library(ottr)"
#RUN R -e "library(here)"
#RUN R -e "library(languageserver)"
# RUN installGithub.r ucbds-infra/ottr@stable

## Install spectratrait and load specific R packages
RUN R -e "devtools::install_github(repo = 'TESTgroup-BNL/spectratrait', dependencies=TRUE)"
RUN R -e "library(here)"

# R Session watcher settings.
# See more details: https://github.com/REditorSupport/vscode-R/wiki/R-Session-watcher
RUN echo 'source(file.path(Sys.getenv("HOME"), ".vscode-R", "init.R"))' >> ${HOME}/.Rprofile

# [Optional] Uncomment this section to install additional OS packages.
# RUN apt-get update \
#     && export DEBIAN_FRONTEND=noninteractive \
#     && apt-get -y install --no-install-recommends <your-package-list-here>

# [Optional] Uncomment this section to install additional R packages.
# RUN install2.r --error --skipinstalled --ncpus -1 <your-package-list-here>

# [Optional] Uncomment this section to install vscode-jupyter dependencies.
#RUN apt-get update \
#    && export DEBIAN_FRONTEND=noninteractive \
#    && apt-get -y install --no-install-recommends libzmq3-dev \
#    && install2.r --error --skipinstalled --ncpus -1 IRkernel \
#    && python3 -m pip --no-cache-dir install jupyter \
#    && R --vanilla -s -e 'IRkernel::installspec(user = FALSE)'

# Expose Rstudio port. Log in with rstudio/rstudio
EXPOSE 8787

# Expose Jupyter port
EXPOSE 8888

CMD ["/init"]