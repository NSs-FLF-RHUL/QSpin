module MFriction

function read_json(file)
    open(file, "r") do f
        return JSON.parse(f)
    end
end

end
