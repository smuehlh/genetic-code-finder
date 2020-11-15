module Statistics
    extend self

    def mean(arr)
        sum = arr.inject(:+)
        len = arr.size
        Statistics.fraction(sum, len)
    end

    def median(arr)
        sorted = arr.sort
        len = sorted.size
        middle = len / 2
        if len % 2 == 0
            Statistics.mean(sorted[middle-1..middle])
        else
            sorted[middle]
        end
    end

    def fraction(num, denom)
        num / denom.to_f
    end

    def percentage(num, total)
        Statistics.fraction(num, total) * 100
    end
end