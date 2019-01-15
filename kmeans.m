K = 10;

X = load('../digit/digit.txt');
Y = load('../digit/labels.txt');
max_iterations = 20;
total_SS_plot = zeros(10, 1);
pairs = zeros(10, 3);

for n=1:10
  for K=1:10
    centroids = initCentroids_random(X, K); %replace by initCentroids_first if required

    for i=1:max_iterations
      prev_centroids = centroids;
      indices = getClosestCentroids(X, centroids);
      centroids = computeCentroids(X, indices, K);
      if prev_centroids == centroids
        %if K == 6
        %    disp("Number of iterations i = " + num2str(i));
        %end
        break
      end
    end

    SS = sum_of_squares(X, indices, centroids);
    total = sum(SS);
    total_SS_plot(K) = total_SS_plot(K) + total;
    
    p1_n = 0;
    p1_d = 0;
    p2_n = 0;
    p2_d = 0;
    m=size(X, 1);
    for xi=1:m
      for xj=xi+1:m
        if Y(xi) == Y(xj)
          p1_d = p1_d + 1;
          if indices(xi) == indices(xj)
            p1_n = p1_n + 1;
          end
        end
        if Y(xi) ~= Y(xj)
          p2_d = p2_d + 1;
          if indices(xi) ~= indices(xj)
            p2_n = p2_n + 1;
          end
        end
      end
    end    
    pairs(K,1) = pairs(K,1) + p1_n/p1_d;
    pairs(K,2) = pairs(K,2) + p2_n/p2_d;
    pairs(K,3) = pairs(K,3) + (p1_n/p1_d + p2_n/p2_d)/2;
  end
end

pairs = pairs / 10; %averaging
pairs = pairs * 100; %percentage
plot(1:10,pairs(:,1))
title("p1");
xlabel("K");
ylabel("p1");
figure()
plot(1:10,pairs(:,2))
title("p2");
xlabel("K");
ylabel("p2");
figure()
plot(1:10,pairs(:,3))
title("p3");
xlabel("K");
ylabel("p3");
figure()
plot(1:10,total_SS_plot/10)
title("Sum of squares");
xlabel("K");
ylabel("Sum of squares");

function centroids = initCentroids_random(X, K)
  centroids = zeros(K,size(X,2)); 
  randidx = randperm(size(X,1));
  centroids = X(randidx(1:K), :);
end

function centroids = initCentroids_first(X, K)
  centroids = zeros(K,size(X,2)); 
  centroids = X(1:K, :);
end

function indices = getClosestCentroids(X, centroids)
  K = size(centroids, 1);
  indices = zeros(size(X,1), 1);
  m = size(X,1);

  for i=1:m
    k = 1;
    min_dist = sum((X(i,:) - centroids(1,:)) .^ 2);
    for j=2:K
        dist = sum((X(i,:) - centroids(j,:)) .^ 2);
        if(dist < min_dist)
          min_dist = dist;
          k = j;
        end
    end
    indices(i) = k;
  end
end

function centroids = computeCentroids(X, idx, K)

  [m n] = size(X);
  centroids = zeros(K, n);
  
  for i=1:K
    xi = X(idx==i,:);
    ck = size(xi,1);
    centroids(i, :) = (1/ck) * sum(xi);
  end
end

function SS = sum_of_squares(X, indices, centroids)
  K = size(centroids, 1);
  m = size(X,1);
  SS = zeros(K,1);
  
  for i=1:m
    SS(indices(i)) = SS(indices(i)) + sum((X(i,:) - centroids(indices(i),:)) .^ 2);
  end
  %SS(m+1) = sum(SS(1:K));
end

