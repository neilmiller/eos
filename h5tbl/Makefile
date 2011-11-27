all : write_idealgas_tbl test_idealh5

test_idealh5 : test_idealh5.f h5table.o ideal_water.o
	gfortran -o test_idealh5 test_idealh5.f h5table.o ideal_water.o -ffree-form

write_idealgas_tbl : write_idealgas_tbl.f ideal_water.o
	gfortran -o write_idealgas_tbl write_idealgas_tbl.f ideal_water.o -ffree-form

ideal_water.o : ideal_water.f
	gfortran -c ideal_water.f -ffree-form

h5table.o : h5table.f
	gfortran -c h5table.f -ffree-form

clean :
	rm ideal_water.o write_idealgas_tbl