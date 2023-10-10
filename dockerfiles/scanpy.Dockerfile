FROM --platform=linux/amd64 gcfntnu/scanpy:1.9.2

MAINTAINER Yuyao Song

CMD ["echo", "Container for harmony and scanorama integration in BENGAL"]

RUN pip3 install harmonypy scanorama click

RUN pip3 install pydantic

ENTRYPOINT ["python"]

