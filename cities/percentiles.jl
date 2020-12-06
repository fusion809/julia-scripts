using Statistics;
using StatsBase;
percentiles = 10:10:90;

function printArrEl(arr)
    """
    Print array elements with proper punctuation between elements

    arr: the array whose elements are to be printed.
    """
    N = length(arr);
    if N == 2
        println(arr[1], " and ", arr[2])
    elseif N == 1
        println(arr[1])
    elseif N != 0
        for i = 1:N-1
            print(arr[i]);
            if i < N-1
                print(", ");
            end
        end
        print(" and ");
        println(arr[N]);
    end
end

function percentilesMsg(cities, data, dataStr)
    """
    Generates a message listing which cities are within each decile for the metric
    in the specified data array whose plain name is in specified string.

    cities:    array of city names as strings.
    data:      numerical data corresponding to each city in cities.
    dataStr:   string stating what data pertains to.
    """
    percentiles = 10:10:90;
    percData = percentile(data, percentiles);

    for i = 1:10
        if i == 1
            decile = data .< percData[i];
            print("1st decile for ", dataStr, " comprises ");
            printArrEl(cities[decile]);
        elseif i == 10
            decile = data .> percData[i-1];
            print("10th decile for ", dataStr, " comprises ");
            printArrEl(cities[decile]);
        elseif i == 2
            decile = percData[i-1] .< data .< percData[i];
            print(i, "nd decile for ", dataStr, " comprises ");
            printArrEl(cities[decile]);
        elseif i == 3
            decile = percData[i-1] .< data .< percData[i];
            print(i, "rd decile for ", dataStr, " comprises ");
            printArrEl(cities[decile]);
        else
            decile = percData[i-1] .< data .< percData[i];
            print(i, "th decile for ", dataStr, " comprises ");
            printArrEl(cities[decile]);
        end
    end

    print("\n");
end


# Aussie cities with universities that have mathematics departments
cities = ["Adelaide"
"Albany"
"Albury"
"Armidale"
"Ballarat"
"Brisbane"
"Canberra"
"Geelong"
"Gold Coast"
"Hobart"
"Lismore"
"Melbourne"
"Newcastle"
"Perth"
"Sydney"
"Toowoomba"
"Townsville"
"Wollongong"];

# Populations in thousands
pops = [1346
34
49
24.5
105
2514
427
263
679
247
28.7
5078
322
2059
5312
137
181
303];

# Minimum recorded temperature
recMin = [-0.4
-0.2
-4
-7
-5.6
2.6
-10
-4.3
2.5
-2.8
-2.9
-2.5
1.8
-0.7
3.7
-1.8
1.1
0.8];

# Annual mean minimum temperature
annMeanMin = [12.3
10.6
9.1
7.5
6.9
16.4
6.8
9
17.3
9
12.9
9.6
14.3
12.9
14.7
12.7
19.8
13.3];

# Mean annual rainfall
annRain = [547.1
798.1
605.7
742.9
620.9
1011.5
583.2
524.3
1259.2
565.6
1162
531.3
1118
733.2
1149.7
699.6
1136
1348.6];

# Maximum record temperature
recMax = [47.7
45.6
46.1
37.1
44.1
41.7
44
47.4
40.5
41.8
42.8
46.8
42.5
44.5
45.8
40.8
44.3
44.1];

# Mean July minimum temperature
julMeanMin = [7.6
7.5
3.2
1.3
3
10.4
0
5.4
12
5.2
6.1
5.5
8.5
7.9
8.9
6.6
13.7
8.3];

# Mean maximum temperature
annMeanMax = [22.5
20.3
22.5
19.6
18
26.6
21.2
19.3
25.4
17.6
26.1
19.9
21.8
24.8
22.8
23.3
29
21.8];

# Mean January maximum temperature
janMeanMax = [29.6
24.8
32.5
26.4
26.2
30.5
30.8
24.5
28.9
22.7
30.4
26.6
25.6
31.2
27
28.5
31.5
25.6];

percentilesMsg(cities, pops, "population");
percentilesMsg(cities, recMin, "minimum record temperature");
percentilesMsg(cities, annMeanMin, "annual average minimum temperature");
percentilesMsg(cities, julMeanMin, "July average minimum temperature");
percentilesMsg(cities, recMax, "maximum recorded temperature");
percentilesMsg(cities, annMeanMax, "annual average maximum temperature");
percentilesMsg(cities, janMeanMax, "January average maximum temperature");
percentilesMsg(cities, annRain, "average annual rainfall");