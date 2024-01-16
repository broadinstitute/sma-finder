TAG1 = weisburd/sma_finder:latest
TAG2 = weisburd/sma_finder:v1.4.1

all: build push update_sha

build:
	cd .. && docker build -f docker/Dockerfile -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)

update_sha:
	docker pull $(TAG1) 2>&1 | grep Digest | cut -c 9- > sha256.txt
	cat sha256.txt && [ ! -z "`cat sha256.txt`" ] && sed -i.bak "s/@sha256:[^\"]*/"@`cat sha256.txt`"/"  ../*.py ../*.wdl
	rm ../*.py.bak ../*.wdl.bak
