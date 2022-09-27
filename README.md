# Project2

To compile the code
```Shell
g++ -I include src/utils.cpp main.cpp -o main
```

And to run it
```Shell
./main > prob5.txt
```

Where `prob5.txt` is a data file that contains all the output, but we will only be using the last part of it. 
See table 1 in the pdf file. 

Running the file also produces the following `.txt` files
```Shell
eig_vec_ana_4.txt  eig_vec_num_4.txt  eig_vec_num_6a.txt  eig_vec_num_6b.txt  prob5.txt  x_axis_6a.txt  x_axis_6b.txt
```

Which is used to the corresponding `plot_eig_vec.py` which can be run as

```Shell
python plot_eig_vec.py
```

To run the test for Problem 3 we can run the command
```Shell
g++ -I include src/utils.cpp test_3.cpp -o test_3
```
