data(fl25)
data(fl25_graph)
data(fl25_enum)
adj = fl25_graph
pop = fl25$pop
plans_10 = fl25_enum$plans[, fl25_enum$pop_dev <= 0.10]

RNGkind("Mersenne-Twister", "Inversion", "Rejection" )
