FROM mcr.microsoft.com/dotnet/runtime:5.0 AS base
WORKDIR /app

FROM mcr.microsoft.com/dotnet/sdk:5.0 AS build
WORKDIR /src
COPY ["Chislaki4/Chislaki4.csproj", "Chislaki4/"]
RUN dotnet restore "Chislaki4/Chislaki4.csproj"
COPY . .
WORKDIR "/src/Chislaki4"
RUN dotnet build "Chislaki4.csproj" -c Release -o /app/build

FROM build AS publish
RUN dotnet publish "Chislaki4.csproj" -c Release -o /app/publish

FROM base AS final
WORKDIR /app
COPY --from=publish /app/publish .
ENTRYPOINT ["dotnet", "Chislaki4.dll"]
