CFLAGS+=-Wall -Wextra

LDFLAGS+=-lm

SAN?=0
SCAN?=0
DEBUG?=0

ifeq ($(SAN), 1)
CFLAGS+=-fsanitize=address
DEBUG=1
endif

ifeq ($(ANALYZER),1)
CFLAGS+=-fanalyzer
DEBUG=1
endif

ifeq ($(DEBUG),1)
CFLAGS+=-g3
else
CFLAGS+=-O3 -DNDEBUG
CFLAGS+=-fno-math-errno # Fine since we don't catch them in any case
CFLAGS+=-ffinite-math-only # Yes, points should be in the domain
CFLAGS+=-fno-signed-zeros # Don't care about that
CFLAGS+=-fno-trapping-math # neither this
LDFLAGS+=-flto
endif

test_collide3_f32: test/collide3_test.c collide3.c
	$(CC) $(CFLAGS) test/collide3_test.c collide3.c $(LDFLAGS) -o test_collide3_f32

test_collide3_f64: test/collide3_test.c collide3.c
	$(CC) $(CFLAGS) -DCOLLIDE3_f64 test/collide3_test.c collide3.c $(LDFLAGS) -o test_collide3_f64
