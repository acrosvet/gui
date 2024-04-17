# layout.jl

using Genie.Renderer.Html

function main_layout(content::String)
  html"""
  <!DOCTYPE html>
  <html lang="en">
  <head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>SARMS-DH/title> <!-- Set your browser tab title here -->
  
  </head>
  <body>
    $content
  </body>
  </html>
  """
end
