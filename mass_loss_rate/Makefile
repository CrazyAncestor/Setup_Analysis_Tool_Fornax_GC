CC = g++
OBJ = main.o

OPENFLAG=-fopenmp


mass_loss_estimater:$(OBJ)
	$(CC) $(OPENFLAG) -o main $(OBJ)  
	rm -f $(OBJ)
	

$(OBJ): %.o: %.cpp
	$(CC) $(OPENFLAG)  -c $<
	
clean:
	rm -f main $(OBJ)
