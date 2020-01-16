function name = makeValidSheetName(name)
    name = name(1:31); %first 31 characters
    name = regexprep(name, '\\|\/|\*|\?|:|\[|\]|\W', '_'); %remove invalid chars - listed & not alphaneumeric
end