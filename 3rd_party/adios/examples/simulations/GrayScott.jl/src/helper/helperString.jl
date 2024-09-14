
export get_type

function get_type(input::String)
    if input == "Float64"
        return Float64
    elseif input == "Float32"
        return Float32
    elseif input == "Float16"
        return Float16
    end

    return nothing
end