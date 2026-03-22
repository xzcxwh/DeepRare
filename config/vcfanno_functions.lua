-- Helper functions for vcfanno to aggregate numeric fields

local function to_number(value)
  if value == nil then
    return nil
  end

  local v = tostring(value)
  if v == "." or v == "" then
    return nil
  end

  return tonumber(v)
end

function max_numeric(values)
  local max_val = nil

  for _, value in ipairs(values) do
    local num = to_number(value)
    if num ~= nil then
      if max_val == nil or num > max_val then
        max_val = num
      end
    end
  end

  if max_val == nil then
    return "."
  end

  return tostring(max_val)
end

function min_numeric(values)
  local min_val = nil

  for _, value in ipairs(values) do
    local num = to_number(value)
    if num ~= nil then
      if min_val == nil or num < min_val then
        min_val = num
      end
    end
  end

  if min_val == nil then
    return "."
  end

  return tostring(min_val)
end
