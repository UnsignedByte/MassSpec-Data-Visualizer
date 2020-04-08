recursiveBinary <- function(n, l) { #Finds binary numbers in string form given length l and n on bits (n<=l)
	# print(paste(n, " ", l))
	if (n==0){
		return(0);
	}else if (n == l){
		return(bitwShiftL(1,l)-1);
	}else{
		return(c(recursiveBinary(n,l-1), bitwShiftL(1,l-1)+recursiveBinary(n-1,l-1)));
	}
}
for(i in 1:4){
	print(recursiveBinary(i,4))
}