setup=mkl

if grep -Fq "FlaringIndex	        0.5" "$setup".par
then

echo "GOOD"
else
echo "BAD"
fi
