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

# Install needed packages and setup non-root user. Use a separate RUN statement to add your own dependencies.
ARG USERNAME=rstudio
ARG USER_UID=1000
ARG USER_GID=$USER_UID
COPY library-scripts/*.sh /tmp/library-scripts/
RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
    && /bin/bash /tmp/library-scripts/common-debian.sh "${INSTALL_ZSH}" "${USERNAME}" "${USER_UID}" "${USER_GID}" "${UPGRADE_PACKAGES}" "true" "true" \
    && usermod -a -G staff ${USERNAME} \
    && apt-get -y install \
        python3-pip \
        libgit2-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libxt-dev \
    && apt-get autoremove -y && apt-get clean -y && rm -rf /var/lib/apt/lists/* /tmp/library-scripts \
    && python3 -m pip --no-cache-dir install radian \
    && install2.r --error --skipinstalled --ncpus -1 \
        devtools (>= 2.3.1) \
        here (>= 0.1) \
        ggplot2 (>= 3.3.2) \
        gridExtra (>= 2.3) \
        remotes (>= 2.2.0) \
        data.table \
        readr (>= 1.3.1) \
        RCurl (>= 1.98-1.2) \
        httr (>= 1.4.2) \
        dplyr (>= 1.0.1) \
        magrittr (>= 2.0.1) \
        plotrix (>= 3.7-8) \
        scales (>= 1.1.1) \
        reshape2 (>= 1.4.4) \
        pls (>= 2.7-2) \
        knitr (>= 1.37) \
        markdown \
        rmarkdown \
        roxygen2 \
        testthat (>= 3.1.2) \
    && rm -rf /tmp/downloaded_packages

## Install summarytools and load some R package off the bat
#RUN R -e "devtools::install_github('https://github.com/dcomtois/summarytools/tree/0-8-9')"
#RUN R -e "library(ottr)"
#RUN R -e "library(here)"
#RUN R -e "library(languageserver)"
# RUN installGithub.r ucbds-infra/ottr@stable

## Install spectratrait and load specific R packages
# or should we leave out most of the packages above and rely on this dependency build to work for the dev env?
RUN R -e "devtools::install_github(repo = 'TESTgroup-BNL/spectratrait', dependencies=TRUE)"
RUN R -e "library(here)"

# VSCode R Debugger dependency. Install the latest release version from GitHub without using GitHub API.
# See https://github.com/microsoft/vscode-dev-containers/issues/1032
RUN export TAG=$(git ls-remote --tags --refs --sort='version:refname' https://github.com/ManuelHentschel/vscDebugger v\* | tail -n 1 | cut --delimiter='/' --fields=3) \
    && Rscript -e "remotes::install_git('https://github.com/ManuelHentschel/vscDebugger.git', ref = '"${TAG}"', dependencies = FALSE)"

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