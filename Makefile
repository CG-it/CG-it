
.PHONY: all install archive update install-sdk update-sdk

VERSION = 1.1
GIT= $(shell which git)

CPP=g++
CPPFLAGS=-fpic -O3 -I${TCLINC}

print-%: ; @echo $*=$($*)

all: libcgmap.so

## Install everything from .rsync-include
install: all
	mkdir -p ${HOME}/.vmdplugins/cgit${VERSION} 2>/dev/null || :
	rsync -av --include-from=.rsync-include --exclude="*" . ${HOME}/.vmdplugins/cgit${VERSION}/.
	cd cgmap; make install "VERSION=$(VERSION)"

## Make archive from latest commit
archive:
	$(GIT) archive --prefix=cgit${VERSION}/ HEAD -o cgit-latest.zip
	cd examples; make archive

## Update list of files to be installed
update:
	rm -rf .rsync-include 2>/dev/null || :
	$(GIT) archive HEAD | tar -t > .rsync-include

## Install the sdk forcefield
install-sdk:
ifneq ($(wildcard ./sdk/.),)
else
	$(GIT) clone git@bitbucket.org:macdercm/sdk.git sdk
endif

# Update sdk forcefield
update-sdk:
ifneq ($(wildcard ./sdk/.),)
	cd ./sdk && $(GIT) pull && make
else
endif

## Make cgMap.so
libcgmap.so: 
	cd cgmap; make "VERSION=$(VERSION)"
