function name = makeValidSheetName(name)
    name = name(1:min(31,end)); %first 31 characters
    name = regexprep(name, '\\|\/|\*|\?|:|\[|\]|\W', '_'); %remove invalid chars - listed & not alphaneumeric
end