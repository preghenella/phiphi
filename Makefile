all: phiphi

phiphi: phiphi.cc
	g++ -o phiphi phiphi.cc \
	-std=c++11 \
	-I${PYTHIA_ROOT}/include -L${PYTHIA_ROOT}/lib -lpythia8 \
	`root-config --cflags --libs` \
	-I${BOOST_ROOT}/include -L${BOOST_ROOT}/lib -lboost_program_options

clean:
	rm -rf phiphi *~

