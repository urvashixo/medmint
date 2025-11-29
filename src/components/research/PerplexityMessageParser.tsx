import React from 'react'

interface PerplexityMessageParserProps {
  content: string
  citations?: string[]
}

export function PerplexityMessageParser({ content, citations }: PerplexityMessageParserProps) {
  // Parse markdown-like formatting from Perplexity responses
  const parseContent = (text: string) => {
    // Split by code blocks first
    const parts = text.split(/(```[\s\S]*?```)/g)
    
    return parts.map((part, index) => {
      if (part.startsWith('```') && part.endsWith('```')) {
        // Code block
        const code = part.slice(3, -3).trim()
        const lines = code.split('\n')
        const language = lines[0] || ''
        const codeContent = lines.slice(1).join('\n') || code
        
        return (
          <div key={index} className="my-4">
            <div className="bg-gray-900 border border-gray-600 rounded-lg overflow-hidden">
              {language && (
                <div className="bg-gray-800 px-3 py-1 text-xs text-gray-400 border-b border-gray-600">
                  {language}
                </div>
              )}
              <pre className="p-4 text-sm text-gray-300 overflow-x-auto">
                <code>{codeContent}</code>
              </pre>
            </div>
          </div>
        )
      } else {
        // Regular text with inline formatting
        return (
          <div key={index} className="whitespace-pre-wrap">
            {formatInlineText(part)}
          </div>
        )
      }
    })
  }

  const formatInlineText = (text: string) => {
    // Handle bold text
    let formatted = text.replace(/\*\*(.*?)\*\*/g, '<strong>$1</strong>')
    
    // Handle italic text
    formatted = formatted.replace(/\*(.*?)\*/g, '<em>$1</em>')
    
    // Handle inline code
    formatted = formatted.replace(/`(.*?)`/g, '<code class="bg-gray-700 px-1 py-0.5 rounded text-sm">$1</code>')
    
    return <span dangerouslySetInnerHTML={{ __html: formatted }} />
  }

  const handleCitationClick = (citation: string) => {
    // Clean the citation string by removing angle brackets and trimming
    let url = citation.trim().replace(/^<+|>+$/g, '')
    
    // If it doesn't start with http/https, add https://
    if (!url.startsWith('http://') && !url.startsWith('https://')) {
      url = 'https://' + url
    }
    
    // Open in new tab
    window.open(url, '_blank', 'noopener,noreferrer')
  }

  return (
    <div className="prose prose-invert max-w-none">
      {parseContent(content)}
      
      {/* Show citations inline if available */}
      {citations && citations.length > 0 && (
        <div className="mt-4 pt-3 border-t border-gray-600">
          <div className="text-xs text-gray-400 mb-2">Sources:</div>
          <div className="space-y-1">
            {citations.slice(0, 3).map((citation, index) => {
              const cleanUrl = citation.trim().replace(/^<+|>+$/g, '')
              const displayUrl = cleanUrl.length > 50 ? cleanUrl.substring(0, 50) + '...' : cleanUrl
              
              return (
                <div key={index} className="text-xs">
                  <button 
                    onClick={() => handleCitationClick(citation)}
                    className="text-blue-400 hover:text-blue-300 underline cursor-pointer"
                  >
                    [{index + 1}] {displayUrl}
                  </button>
                </div>
              )
            })}
            {citations.length > 3 && (
              <div className="text-xs text-gray-500">
                +{citations.length - 3} more sources
              </div>
            )}
          </div>
        </div>
      )}
    </div>
  )
}