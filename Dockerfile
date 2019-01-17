FROM kbase/sdkbase2:python
MAINTAINER KBase Developer
ENV SAMTOOLS_VER 1.9

# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

# RUN apt-get update

# Install Bowtie2.2.9 (Bowtie2 in KBase 2019-01-16 is version 2.3.2)

# Install Bowtie2.2.9 (Bowtie2 in KBase 2019-01-16 is version 2.3.2)
RUN VERSION='2.2.9' \
    && mkdir bowtie2-bin \
    && curl -L -o bowtie2-${VERSION}.zip "https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/bowtie2-2.2.9-linux-x86_64.zip/download" \
    && unzip -q bowtie2-${VERSION}.zip \
    && cd bowtie2-${VERSION} \
    && cp bowtie2 bowtie2-align-l bowtie2-align-s bowtie2-build bowtie2-build-l bowtie2-build-s \
          bowtie2-inspect bowtie2-inspect-l bowtie2-inspect-s ../bowtie2-bin \
    && cd .. \
    && rm -rf bowtie2-${VERSION}*

# Install samtools
RUN conda install -c bioconda samtools

ENV PATH $PATH:/bowtie2-bin


# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
