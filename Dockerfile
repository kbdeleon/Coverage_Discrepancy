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
RUN VERSION='2.2.9' \
    && mkdir bowtie2-bin \
    && wget --no-verbose "http://sourceforge.net/projects/bowtie-bio/files/bowtie2/${VERSION}/bowtie2-${VERSION}-source.zip" \
    && unzip -q bowtie2-${VERSION}-source.zip \
    && cd bowtie2-${VERSION} \
    && make NO_TBB=1 \
    && cp bowtie2 bowtie2-align-l bowtie2-align-s bowtie2-build bowtie2-build-l bowtie2-build-s \
          bowtie2-inspect bowtie2-inspect-l bowtie2-inspect-s ../bowtie2-bin \
    && cd .. \
    && rm -rf bowtie2-${VERSION}*

# Install samtools
RUN cd /opt \
    && wget https://github.com/samtools/samtools/releases/download/$SAMTOOLS_VER/samtools-$SAMTOOLS_VER.tar.bz2 \
    && tar xvjf samtools-$SAMTOOLS_VER.tar.bz2 \
    && rm -f samtools-$SAMTOOLS_VER.tar.bz2 \
    && cd samtools-$SAMTOOLS_VER \
    && ./configure \
    && make \
    && make install

ENV PATH $PATH:/opt/samtools-$SAMTOOLS_VER



# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
