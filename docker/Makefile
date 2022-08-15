TAG = weisburd/sma_finder:latest

all: build push update_sha

build:
	cd .. && docker build -f docker/Dockerfile -t $(TAG) .

push:
	docker push $(TAG)

update_sha:
	docker pull $(TAG) 2>&1 | grep Digest | cut -c 9- > sha256.txt
	cat sha256.txt && [ ! -z "`cat sha256.txt`" ] && sed -i.bak "s/@sha256:[^\"]*/"@`cat sha256.txt`"/"  ../*.py
