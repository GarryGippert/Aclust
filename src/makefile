BIN = $(HOME)/bin/`arch`
#BIN = /work3/garryg/bioC/bin/`arch`
VERSION = `git log -n 1 --pretty="format:%h %ci"`
LFLAGS = -lm
EXE = aclust bclust

.SUFFIXES:
.SUFFIXES: .c .o
 
.c:
	$(CC) $? -o $@ $(LIB) $(LFLAGS)
 
.c.o:   $<
	$(CC) -c $<

all: $(EXE)

save:
	echo "#define VERSION \"$(VERSION)\"" > version.txt
	git commit aclust.c bclust.c version.txt makefile
	git push

clean:
	/bin/rm -rf *.o *~ core* $(EXE)

install: $(EXE)
	@-files="$(EXE)"; \
	for f in $$files; do \
	  target=$(BIN)/`basename $$f`; \
	  install -m 775 $$f $$target; \
	  /bin/rm -f $$f; \
	  echo $$f installed at $$target and removed locally; \
	done
