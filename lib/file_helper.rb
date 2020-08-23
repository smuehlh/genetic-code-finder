module FileHelper
    extend self

    def file_exist_or_die(path)
        if ! FileTest.file?(path) then
            abort "File #{path} does not exist."
        end
    end
end