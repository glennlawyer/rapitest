FROM trestletech/plumber

MAINTAINER Glenn Lawyer <dr.g.lawyer@gmail.com>

# Create workdir and user
RUN mkdir /app && groupadd -r user && useradd -r -g user user && mkdir -p /home/user/r_libs && chown -R user:user /home/user
ENV R_LIBS /home/user/r_libs
WORKDIR /app
# Check your privileges, yo!!
USER user

# install R dependencies
RUN R -e 'install.packages(c("rdist", "data.table", "Matrix"), repos = "http://cran.us.r-project.org")'

## copy the app into the container
COPY clusterMoods.R /app/

## run it (via plumber)
CMD ["/app/clusterMoods.R"]
