function k = resample2(w)

[w,I] = sort(w,'descend');
N = length(w);
N = N/2;
Q = cumsum(w);
T = zeros(1,N+1);
index = zeros(1,N);

T = linspace(0,1-1/N,N) + rand()/N;
T(N+1) = 1;

i = 1;
j = 1;

while(i<=N)
    if (T(i)<Q(j))
        index(i) = j;
        i = i + 1;
    else
        j = j + 1;
    end
end

k = I(index);