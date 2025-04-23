# 0 31 32 63
method=$1
batch=$2

for i in `seq 0 ${batch}`; do
	nohup python scipt/function_${method}.py ${batch} ${i} &
	#echo ${method}_${i}
done

